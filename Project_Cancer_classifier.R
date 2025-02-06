#Related packages
library(caret)
library(pROC)
library(DALEX)
library(dplyr)

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("TCGAbiolinks")
BiocManager::install("SummarizedExperiment")

library(TCGAbiolinks)
library(SummarizedExperiment)

##Query TCGA for methylation data

query <- GDCquery(
  project = "TCGA-BRCA",
  data.category = "DNA Methylation",
  data.type = "Methylation Beta Value",
  platform = "Illumina Human Methylation 27")

brca_res = getResults(query)
summary(factor(brca_res$sample_type))

#filter sample_type with less number
query <- GDCquery(
  project = "TCGA-BRCA",
  data.category = "DNA Methylation",
  data.type = "Methylation Beta Value",
  platform = "Illumina Human Methylation 27",
  sample.type = c("Primary Tumor", "Solid Tissue Normal"))

##Download the data
GDCdownload(query)

##Prepare and load the data
methylation_data <- GDCprepare(query)

##Get patient metadata
clinical_metadata <- as.data.frame(colData(methylation_data))
clinical_selected <- clinical_metadata %>% select(sample_type)
methylation_beta <- as.data.frame(assay(methylation_data))
#Remove rows with NA values
methylation_beta_clean <- na.omit(methylation_beta)

#Preprocess the data
#Select top variable genes
all.trans <- data.frame(t(methylation_beta_clean))
SDs = apply(all.trans,2, sd)
topPreds = order(SDs, decreasing = TRUE)[1:1000]
all.trans=all.trans[,topPreds]
#Add metadata information
all.trans <- merge(all.trans,clinical_selected, by = "row.names")
rownames(all.trans) <- all.trans$Row.names
all.trans <- all.trans[,-1]
#Remove Near Zero variation
all.zero <- preProcess(all.trans, method = 'nzv',uniqueCut = 15)
all.trans <- predict(all.zero, all.trans)
#Center (Normalization)
all.center <- preProcess(all.trans, method = 'center')
all.trans <- predict(all.center, all.trans)
#Remove highly correlated
all.corr <- preProcess(all.trans, method = 'corr', cutoff=0.5)
all.trans <- predict(all.corr, all.trans)

##Split the Data into training and Test Sets
intrain <- createDataPartition(y=all.trans$sample_type,p=0.7)[[1]]
train.log <- all.trans[intrain,]
test.log <- all.trans[-intrain,]
#Train a k-Nearest (kNN) Neighbors (kNN) Model
ctrl.lgg <- trainControl(method = 'cv',number=5)
knn.lgg <- train(sample_type~.,data=train.log, method='knn',trControl=ctrl.lgg,tuneGrid=data.frame(k=1:20))
#Predict
trainPred <- predict(knn.lgg,newdata =train.log)
testPred <- predict(knn.lgg,newdata =test.log)
#Interpretation
#Confusion matrix
confusionMatrix(trainPred,factor(train.log$sample_type))
confusionMatrix(testPred,factor(test.log$sample_type))