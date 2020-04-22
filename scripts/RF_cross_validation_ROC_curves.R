library(caret)
library(randomForest)
library(ROCR)
source("/hps/nobackup/ma/natalja/data/mQTL_data/scripts/load_data.R")

annotated_metabolites <-  c("6.30_602.2977m/z","6.13_685.2611m/z", "6.13_600.2850m/z", "6.13_600.4780m/z", "6.13_668.2724m/z", "3.80_258.1087m/z", "4.55_190.2275n", "3.80_126.0658m/z", "8.01_639.3614n", "2.62_125.0837n", "6.13_736.2562m/z", "2.62_92.1133n", "2.03_112.0501m/z", "X5.4085", "3.74_249.0860n", "X1.2051", "X2.873", "X1.7312", "X5.2831", "X5.276", "X3.0968", "X3.0962", "X1.1952", "X3.0979", "X3.0973")

annotated_metabolites2 <-  c("6.30_602.2977m/z","6.13_685.2611m/z", "6.13_600.2850m/z", "6.13_600.4780m/z", "6.13_668.2724m/z", "3.80_258.1087m/z", "4.55_190.2275n", "3.80_126.0658m/z", "8.01_639.3614n", "2.62_125.0837n", "6.13_736.2562m/z", "2.62_92.1133n", "2.03_112.0501m/z", "X5.4085", "3.74_249.0860n", "X1.2051", "X1.1952")

## Select 24 annotated metabolites
# 1) normalised data
met_norm<- expr[rownames(expr) %in% annotated_metabolites2,]
# 2) original data
met_orig<- expr_orig[rownames(expr_orig) %in% annotated_metabolites2,]
# 3) original QN data
met_orig_qn<- expr_orig_qn[rownames(expr_orig_qn) %in% annotated_metabolites2,]

dim(met_orig)
dim(met_norm)
dim(met_orig_qn)
# both dim should be 24x471

cvrt = read.table("/hps/nobackup/ma/natalja/data/mQTL_data/covariates/Covariates_all_vertical_text2.txt",sep="\t",header=T,row.names=1)
cvrt <- cvrt[rownames(cvrt) %in% colnames(expr),]
cvrt <- cvrt[colnames(expr),]

############## Load feature set: metabolites + covariates ######################
LOAD_FS <- function(met){
  df <- as.data.frame(t(met))
  df$Gender <- cvrt$Gender
  df$Centre <- cvrt$Centre
  df$Age <- cvrt$Age
  df$Diagnosis <- cvrt$Diagnosis
  data<- df[cvrt$Diagnosis %in% c("CTL","ADC"),]
  data$Diagnosis<- as.factor(data$Diagnosis)
  names(data) <- make.names(names(data)) 
  data$Diagnosis <- relevel(data$Diagnosis,ref="CTL")
  data$Diagnosis <- as.factor(as.character(data$Diagnosis))
  #data<-data[sample(nrow(data)),]
  return(data)
}  

LOAD_FS_MCI <- function(met){
  df <- as.data.frame(t(met))
  df$Gender <- cvrt$Gender
  df$Centre <- cvrt$Centre
  df$Age <- cvrt$Age
  df$Diagnosis <- cvrt$Diagnosis
  data<- df[cvrt$Diagnosis %in% c("cMCI","sMCI"),]
  data$Diagnosis<- as.factor(data$Diagnosis)
  names(data) <- make.names(names(data)) 
  data$Diagnosis <- relevel(data$Diagnosis,ref="sMCI")
  data$Diagnosis <- as.factor(as.character(data$Diagnosis))
  #data<-data[sample(nrow(data)),]
  return(data)
}  
################################################################################
data_norm <- LOAD_FS(met_norm)
data_orig <- LOAD_FS(met_orig)
data_orig_qn <- LOAD_FS(met_orig_qn)

# Split data 70/30
split=0.70
trainIndex <- createDataPartition(data_norm$Diagnosis, p=split, list=FALSE)

trainData_norm <- data_norm[trainIndex,]
testData_norm <- data_norm[-trainIndex,]

trainData_orig <- data_orig[ trainIndex,]
testData_orig <- data_orig[-trainIndex,]

trainData_orig_qn <- data_orig_qn[ trainIndex,]
testData_orig_qn <- data_orig_qn[-trainIndex,]

# 1) Normalised data
RF_train=randomForest(Diagnosis ~ ., data = trainData_norm, importance = TRUE, proximity=TRUE)
# OOB 7.21%
model_norm <- RF_train


model_lr_norm <- glm(Diagnosis ~ .,family=binomial(link='logit'),data=trainData_norm)

# 2) Original data
# RF
RF_train=randomForest(Diagnosis ~ ., data = trainData_orig, importance = TRUE, proximity=TRUE)
# OOB 33.17%
model_orig <- RF_train
model_lr_orig <- glm(Diagnosis ~ .,family=binomial(link='logit'),data=trainData_orig)

# 3) Original QN data
RF_train=randomForest(Diagnosis ~ ., data = trainData_orig_qn, importance = TRUE, proximity=TRUE)
# OOB 31.73%
model_orig_qn <- RF_train
model_lr_orig_qn <- glm(Diagnosis ~ .,family=binomial(link='logit'),data=trainData_orig_qn)

############### ROC curves RF ##################################################
# 1) Normalised data
RF_predict = predict(model_norm, type = "prob", newdata = testData_norm)
RF_prediction = prediction(RF_predict[,2], as.numeric(testData_norm$Diagnosis)-1)
RF_performance = performance(RF_prediction, "tpr", "fpr")
ROC_norm <- roc(as.numeric(testData_norm$Diagnosis)-1,RF_predict[,2])

perf_AUC=performance(RF_prediction,"auc") #Calculate the AUC value
perf_AUC@y.values[[1]]
#  0.8413462

# 2) Original data
RF_predict = predict(model_orig, type = "prob", newdata = testData_orig)
RF_prediction = prediction(RF_predict[,2], as.numeric(testData_orig$Diagnosis)-1)

RF_performance = performance(RF_prediction, "tpr", "fpr")
ROC_orig <- roc(as.numeric(testData_orig$Diagnosis)-1,RF_predict[,2])

perf_AUC=performance(RF_prediction,"auc") #Calculate the AUC value
perf_AUC@y.values[[1]]

png("ROC_24_metabolites_RF2.png", width = 4, height = 4, units = 'in', res=300)
plot(smooth(ROC_norm),col="blue")
plot(smooth(ROC_orig),col="red", add=TRUE)
legend("bottomright", legend=c("Normalised", "Original"),
       col=c("blue", "red", "orange"), lty=c(1,1,1), cex=0.8,
       box.lty=0)
dev.off()

# # 3) Original QN data
# data <- LOAD_FS(met_orig_qn)
# #split=0.70
# #trainIndex <- createDataPartition(data$Diagnosis, p=split, list=FALSE)
# trainData <- data[ trainIndex,]
# testData <- data[-trainIndex,]
# 
# RF_predict = predict(model_orig_qn, type = "prob", newdata = testData)
# RF_prediction = prediction(RF_predict[,2], as.numeric(testData$Diagnosis)-1)
# RF_performance = performance(RF_prediction, "tpr", "fpr")
# ROC <- roc(as.numeric(testData$Diagnosis)-1,RF_predict[,2])
# 
# perf_AUC=performance(RF_prediction,"auc") #Calculate the AUC value
# perf_AUC@y.values[[1]]
# #0.968956
# plot(smooth(ROC),col="orange", add=TRUE)



############### ROC curves Logistic Regression ##################################################
# 1) Normalised data
RF_predict = predict(model_lr_norm, type = "response", newdata = testData_norm)
RF_prediction = prediction(RF_predict, as.numeric(testData_norm$Diagnosis)-1)
RF_performance = performance(RF_prediction, "tpr", "fpr")
ROC_norm_lr <- roc(as.numeric(testData_norm$Diagnosis)-1,RF_predict)

perf_AUC=performance(RF_prediction,"auc") #Calculate the AUC value
perf_AUC@y.values[[1]]

# 2) Original data
RF_predict = predict(model_lr_orig, type = "response", newdata = testData_orig)
RF_prediction = prediction(RF_predict, as.numeric(testData_orig$Diagnosis)-1)
RF_performance = performance(RF_prediction, "tpr", "fpr")
ROC_orig_lr <- roc(as.numeric(testData_orig$Diagnosis)-1,RF_predict)

perf_AUC=performance(RF_prediction,"auc") #Calculate the AUC value
perf_AUC@y.values[[1]]


png("ROC_24_metabolites_LR2.png", width = 4, height = 4, units = 'in', res=300)
plot(smooth(ROC_norm_lr),col="blue")
plot(smooth(ROC_orig_lr),col="red", add=TRUE)
legend("bottomright", legend=c("Normalised", "Original"),
       col=c("blue", "red", "orange"), lty=c(1,1,1), cex=0.8,
       box.lty=0)
dev.off()

# # 3) Original QN data
# data <- LOAD_FS(met_orig_qn)
# #split=0.70
# #trainIndex <- createDataPartition(data$Diagnosis, p=split, list=FALSE)
# trainData <- data[ trainIndex,]
# testData <- data[-trainIndex,]
# 
# RF_predict = predict(model_lr_orig_qn, type = "response", newdata = testData)
# RF_prediction = prediction(RF_predict, as.numeric(testData$Diagnosis)-1)
# RF_performance = performance(RF_prediction, "tpr", "fpr")
# ROC <- roc(as.numeric(testData$Diagnosis)-1,RF_predict)
# 
# perf_AUC=performance(RF_prediction,"auc") #Calculate the AUC value
# perf_AUC@y.values[[1]]
# #0.968956
# plot(smooth(ROC),col="orange", add=TRUE)



############### ROC curves RF MCIs ##################################################
# 1) Normalised data
data <- LOAD_FS_MCI(met_norm)

RF_predict = predict(model_norm, type = "prob", newdata = data)
RF_prediction = prediction(RF_predict[,2], as.numeric(data$Diagnosis)-1)

RF_performance = performance(RF_prediction, "tpr", "fpr")
ROC_norm_MCI <- roc(as.numeric(data$Diagnosis)-1,RF_predict[,2])

perf_AUC=performance(RF_prediction,"auc") #Calculate the AUC value
perf_AUC@y.values[[1]]


# 2) Original data
data <- LOAD_FS_MCI(met_orig)

RF_predict = predict(model_orig, type = "prob", newdata = data)
RF_prediction = prediction(RF_predict[,2], as.numeric(data$Diagnosis)-1)
RF_performance = performance(RF_prediction, "tpr", "fpr")
ROC_orig_MCI <- roc(as.numeric(data$Diagnosis)-1,RF_predict[,2])

perf_AUC=performance(RF_prediction,"auc") #Calculate the AUC value
perf_AUC@y.values[[1]]

# 0.9865385
png("ROC_24_metabolites_RF_MCIs.png", width = 4, height = 4, units = 'in', res=300)
plot(smooth(ROC_norm_MCI),col="blue")
plot(smooth(ROC_orig_MCI),col="red", add=TRUE)
legend("bottomright", legend=c("Normalised", "Original"),
       col=c("blue", "red", "orange"), lty=c(1,1,1), cex=0.8,
       box.lty=0)
dev.off()
# # 3) Original QN data
# data <- LOAD_FS_MCI(met_orig_qn)
# # #split=0.70
# # #trainIndex <- createDataPartition(data$Diagnosis, p=split, list=FALSE)
# # trainData <- data[ trainIndex,]
# # testData <- data[-trainIndex,]
# # 
# RF_predict = predict(model_orig_qn, type = "prob", newdata = data)
# RF_prediction = prediction(RF_predict[,2], as.numeric(data$Diagnosis)-1)
# RF_performance = performance(RF_prediction, "tpr", "fpr")
# ROC <- roc(as.numeric(data$Diagnosis)-1,RF_predict[,2])
# # 
# perf_AUC=performance(RF_prediction,"auc") #Calculate the AUC value
# perf_AUC@y.values[[1]]
# #0.968956
# plot(smooth(ROC),col="orange", add=TRUE)



############### ROC curves MCI build ##################################################
# 1) Normalised data
data <- LOAD_FS_MCI(met_norm)
split=0.70
trainIndex <- createDataPartition(data$Diagnosis, p=split, list=FALSE)
trainData <- data[ trainIndex,]
testData <- data[-trainIndex,]

RF_train=randomForest(Diagnosis ~ ., data = trainData, importance = TRUE, proximity=TRUE)

RF_predict = predict(RF_train, type = "prob", newdata = testData)
RF_prediction = prediction(RF_predict[,2], as.numeric(testData$Diagnosis)-1)
RF_performance = performance(RF_prediction, "tpr", "fpr")
ROC <- roc(as.numeric(testData$Diagnosis)-1,RF_predict[,2])

perf_AUC=performance(RF_prediction,"auc") #Calculate the AUC value
perf_AUC@y.values[[1]]
# 0.9865385
png("ROC_24_metabolites_RF2.png", width = 4, height = 4, units = 'in', res=300)
plot(smooth(ROC),col="blue")

# 2) Original data
data <- LOAD_FS_MCI(met_orig)
#split=0.70
#trainIndex <- createDataPartition(data$Diagnosis, p=split, list=FALSE)
trainData <- data[ trainIndex,]
testData <- data[-trainIndex,]

RF_train=randomForest(Diagnosis ~ ., data = trainData, importance = TRUE, proximity=TRUE)

RF_predict = predict(RF_train, type = "prob", newdata = testData)
RF_prediction = prediction(RF_predict[,2], as.numeric(testData$Diagnosis)-1)

RF_performance = performance(RF_prediction, "tpr", "fpr")
ROC <- roc(as.numeric(testData$Diagnosis)-1,RF_predict[,2])

perf_AUC=performance(RF_prediction,"auc") #Calculate the AUC value
perf_AUC@y.values[[1]]
# 0.9513963
plot(smooth(ROC),col="red", add=TRUE)

# # 3) Original QN data
# data <- LOAD_FS(met_orig_qn)
# #split=0.70
# #trainIndex <- createDataPartition(data$Diagnosis, p=split, list=FALSE)
# trainData <- data[ trainIndex,]
# testData <- data[-trainIndex,]
# 
# RF_predict = predict(model_orig_qn, type = "prob", newdata = testData)
# RF_prediction = prediction(RF_predict[,2], as.numeric(testData$Diagnosis)-1)
# RF_performance = performance(RF_prediction, "tpr", "fpr")
# ROC <- roc(as.numeric(testData$Diagnosis)-1,RF_predict[,2])
# 
# perf_AUC=performance(RF_prediction,"auc") #Calculate the AUC value
# perf_AUC@y.values[[1]]
# #0.968956
# plot(smooth(ROC),col="orange", add=TRUE)

legend("bottomright", legend=c("Normalised", "Original"),
       col=c("blue", "red", "orange"), lty=c(1,1,1), cex=0.8,
       box.lty=0)
dev.off()

#######################################Miscellaneous############################
split=0.70
trainIndex <- createDataPartition(data$Diagnosis, p=split, list=FALSE)
trainData <- data[ trainIndex,]
testData <- data[-trainIndex,]

#Perform 10 fold cross validation
#Create 10 equally size folds
folds <- cut(seq(1,nrow(data)),breaks=10,labels=FALSE)
res_set <- data.frame(nr = c(1:10))

df.x <- matrix(NA,nrow=30,ncol=10)
df.y <- matrix(NA,nrow=30,ncol=10)
df.alpha <- matrix(NA,nrow=30,ncol=10)
df.AUC <-c(0)

i=1
testIndexes <- which(folds==i,arr.ind=TRUE)
testData <- data[testIndexes, ]
trainData <- data[-testIndexes, ]
# ntree=680, mtry=90, 
RF_train=randomForest(Diagnosis ~ ., data = trainData, importance = TRUE, proximity=TRUE)
RF_predict = predict(RF_train, type = "prob", newdata = testData)
RF_prediction = prediction(RF_predict[,2], as.numeric(testData$Diagnosis)-1)

RF_performance = performance(RF_prediction, "tpr", "fpr")
ROC <- roc(as.numeric(testData$Diagnosis)-1,RF_predict[,2])

perf_AUC=performance(RF_prediction,"auc") #Calculate the AUC value
df.AUC[i] <- perf_AUC@y.values[[1]]

png("ROC_24_metabolites.png", width = 4, height = 4, units = 'in', res=300)
#par(mar=c(5.1,4.1,6,2.1), xpd=TRUE)
#plot(RF_performance)
plot(smooth(ROC),col="blue")

for(i in 2:10){
  #Segement your data by fold using the which() function 
  testIndexes <- which(folds==i,arr.ind=TRUE)
  testData <- data[testIndexes, ]
  trainData <- data[-testIndexes, ]
  
  split=0.90
  trainIndex <- createDataPartition(data$Diagnosis, p=split, list=FALSE)
  trainData <- data[ trainIndex,]
  testData <- data[-trainIndex,]
  # ntree=680, mtry=90, 
  RF_train=randomForest(Diagnosis ~ ., data = trainData, importance = TRUE, proximity=TRUE)
  RF_predict = predict(RF_train, type = "prob", newdata = testData,predict.all=FALSE)
  RF_prediction = prediction(RF_predict[,2], as.numeric(testData$Diagnosis)-1)
  
  RF_performance = performance(RF_prediction, "tpr", "fpr")
  ROC <- roc(as.numeric(testData$Diagnosis)-1,RF_predict[,2])
  
  
  df.x[1:length(as.vector(RF_performance@x.values)[[1]]),i] <- as.vector(RF_performance@x.values)[[1]]
  df.y[1:length(as.vector(RF_performance@y.values)[[1]]),i] <- as.vector(RF_performance@y.values)[[1]]
  df.alpha[1:length(as.vector(RF_performance@alpha.values)[[1]]),i] <- as.vector(RF_performance@alpha.values)[[1]]
  
  perf_AUC=performance(RF_prediction,"auc") #Calculate the AUC value
  df.AUC[i] <- perf_AUC@y.values[[1]]
  #res_set <- cbind(res_set,perf_AUC@y.values[[1]])
  
  #plot(RF_performance, add = TRUE)
  plot(smooth(ROC),col="blue", add=TRUE)
}
#abline(a=0, b= 1)
dev.off()

