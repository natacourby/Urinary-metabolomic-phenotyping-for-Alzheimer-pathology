library(caret)
library(randomForest)
library(ROCR)
library(pROC)
source("./scripts/feature_selection.R")

annotated_metabolites <- read.csv("./results/annotated_metabolites.txt",sep="\t", header=T)
covariates <- c("Gender","Centre","Age","Diagnosis")

data <- feature_selection_data_preparation("D","III")
data_MCI <- feature_selection_data_preparation("D","IV")

## Select intensities of annotated metabolites
data_annotated<- data[,colnames(data) %in% c(as.character(annotated_metabolites$Metabolic_features),covariates)]
data_annotated_MCI<- data_MCI[,colnames(data_MCI) %in% c(as.character(annotated_metabolites$Metabolic_features),covariates)]

################################################################################
# Split data 70/30
split=0.70
trainIndex <- createDataPartition(data_annotated$Diagnosis, p=split, list=FALSE)

trainData <- data_annotated[trainIndex,]
testData <- data_annotated[-trainIndex,]

# Train
RF_train=randomForest(Diagnosis ~ ., data = trainData, importance = TRUE, proximity=TRUE)
# OOB 3.85%
model_rf_annotated <- RF_train

model_lr_annotated <- glm(Diagnosis ~ .,family=binomial(link='logit'),data=trainData)

############### ROC curves RF ##################################################
RF_predict = predict(model_rf_annotated, type = "prob", newdata = testData)
RF_prediction = prediction(RF_predict[,2], as.numeric(testData$Diagnosis)-1)
RF_performance = performance(RF_prediction, "tpr", "fpr")
ROC_res <- roc(as.numeric(testData$Diagnosis)-1,RF_predict[,2])
perf_AUC=performance(RF_prediction,"auc") #Calculate the AUC value
perf_AUC@y.values[[1]]
# 0.9983974

# MCI cases
RF_predict_MCI = predict(model_rf_annotated, type = "prob", newdata = data_annotated_MCI)
RF_prediction_MCI = prediction(RF_predict_MCI[,2], as.numeric(data_annotated_MCI$Diagnosis)-1)
RF_performance_MCI = performance(RF_prediction_MCI, "tpr", "fpr")
ROC_res_MCI <- roc(as.numeric(data_annotated_MCI$Diagnosis)-1,RF_predict_MCI[,2])
perf_AUC_MCI=performance(RF_prediction_MCI,"auc") #Calculate the AUC value
perf_AUC_MCI@y.values[[1]]
# 0.8691468

png("./images/ROC_RF_annotated_metabolites.png", width = 4, height = 4, units = 'in', res=300)
plot(smooth(ROC_res),col="blue")
plot(smooth(ROC_res_MCI),col="red", add=TRUE)
legend("bottomright", legend=c("AD/CTL", "cMCI/sMCI"),
       col=c("blue", "red", "orange"), lty=c(1,1,1), cex=0.8,
       box.lty=0)
dev.off()
############### ROC curves Logistic Regression ##################################################
RF_predict = predict(model_lr_annotated, type = "response", newdata = testData)
RF_prediction = prediction(RF_predict, as.numeric(testData$Diagnosis)-1)
RF_performance = performance(RF_prediction, "tpr", "fpr")
ROC_res_lr <- roc(as.numeric(testData$Diagnosis)-1,RF_predict)

perf_AUC=performance(RF_prediction,"auc") #Calculate the AUC value
perf_AUC@y.values[[1]]
#0.9941239

# MCI cases
RF_predict_MCI = predict(model_lr_annotated, type = "response", newdata = data_annotated_MCI)
RF_prediction_MCI = prediction(RF_predict_MCI, as.numeric(data_annotated_MCI$Diagnosis)-1)
RF_performance_MCI = performance(RF_prediction_MCI, "tpr", "fpr")
ROC_res_lr_MCI <- roc(as.numeric(data_annotated_MCI$Diagnosis)-1,RF_predict_MCI)

perf_AUC_MCI=performance(RF_prediction_MCI,"auc") #Calculate the AUC value
perf_AUC_MCI@y.values[[1]]

png("./images/ROC_LR_annotated_metabolites.png", width = 4, height = 4, units = 'in', res=300)
plot(smooth(ROC_res_lr),col="blue")
plot(smooth(ROC_res_lr_MCI),col="red", add=TRUE)
legend("bottomright", legend=c("AD/CTL", "cMCI/sMCI"),
       col=c("blue", "red", "orange"), lty=c(1,1,1), cex=0.8,
       box.lty=0)
dev.off()

##### Remove drugs from annotated metabolites

annotated_metabolites <- annotated_metabolites[!annotated_metabolites$Drug,]
## Select intensities of annotated metabolites
data_annotated<- data[,colnames(data) %in% c(as.character(annotated_metabolites$Metabolic_features),covariates)]
data_annotated_MCI<- data_MCI[,colnames(data_MCI) %in% c(as.character(annotated_metabolites$Metabolic_features),covariates)]

################################################################################
# Split data 70/30
split=0.70
trainIndex <- createDataPartition(data_annotated$Diagnosis, p=split, list=FALSE)

trainData <- data_annotated[trainIndex,]
testData <- data_annotated[-trainIndex,]

# Train
RF_train=randomForest(Diagnosis ~ ., data = trainData, importance = TRUE, proximity=TRUE)
# OOB 3.85%
model_rf_annotated <- RF_train

model_lr_annotated <- glm(Diagnosis ~ .,family=binomial(link='logit'),data=trainData)

############### ROC curves RF ##################################################
RF_predict = predict(model_rf_annotated, type = "prob", newdata = testData)
RF_prediction = prediction(RF_predict[,2], as.numeric(testData$Diagnosis)-1)
RF_performance = performance(RF_prediction, "tpr", "fpr")
ROC_res <- roc(as.numeric(testData$Diagnosis)-1,RF_predict[,2])
perf_AUC=performance(RF_prediction,"auc") #Calculate the AUC value
perf_AUC@y.values[[1]]
# 0.920406

# MCI cases
RF_predict_MCI = predict(model_rf_annotated, type = "prob", newdata = data_annotated_MCI)
RF_prediction_MCI = prediction(RF_predict_MCI[,2], as.numeric(data_annotated_MCI$Diagnosis)-1)
RF_performance_MCI = performance(RF_prediction_MCI, "tpr", "fpr")
ROC_res_MCI <- roc(as.numeric(data_annotated_MCI$Diagnosis)-1,RF_predict_MCI[,2])
perf_AUC_MCI=performance(RF_prediction_MCI,"auc") #Calculate the AUC value
perf_AUC_MCI@y.values[[1]]
# 0.7411706

png("./images/ROC_RF_annotated_metabolites_drugs_removed.png", width = 4, height = 4, units = 'in', res=300)
plot(smooth(ROC_res),col="blue")
plot(smooth(ROC_res_MCI),col="red", add=TRUE)
legend("bottomright", legend=c("AD/CTL", "cMCI/sMCI"),
       col=c("blue", "red", "orange"), lty=c(1,1,1), cex=0.8,
       box.lty=0)
dev.off()
############### ROC curves Logistic Regression ##################################################
RF_predict = predict(model_lr_annotated, type = "response", newdata = testData)
RF_prediction = prediction(RF_predict, as.numeric(testData$Diagnosis)-1)
RF_performance = performance(RF_prediction, "tpr", "fpr")
ROC_res_lr <- roc(as.numeric(testData$Diagnosis)-1,RF_predict)

perf_AUC=performance(RF_prediction,"auc") #Calculate the AUC value
perf_AUC@y.values[[1]]
#0.8191774

# MCI cases
RF_predict_MCI = predict(model_lr_annotated, type = "response", newdata = data_annotated_MCI)
RF_prediction_MCI = prediction(RF_predict_MCI, as.numeric(data_annotated_MCI$Diagnosis)-1)
RF_performance_MCI = performance(RF_prediction_MCI, "tpr", "fpr")
ROC_res_lr_MCI <- roc(as.numeric(data_annotated_MCI$Diagnosis)-1,RF_predict_MCI)

perf_AUC_MCI=performance(RF_prediction_MCI,"auc") #Calculate the AUC value
perf_AUC_MCI@y.values[[1]]
#0.6974206

png("./images/ROC_LR_annotated_metabolites_drugs_removed.png", width = 4, height = 4, units = 'in', res=300)
plot(smooth(ROC_res_lr),col="blue")
plot(smooth(ROC_res_lr_MCI),col="red", add=TRUE)
legend("bottomright", legend=c("AD/CTL", "cMCI/sMCI"),
       col=c("blue", "red", "orange"), lty=c(1,1,1), cex=0.8,
       box.lty=0)
dev.off()

