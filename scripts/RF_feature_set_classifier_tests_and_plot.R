################################################################################
# RF first run for all feature set and classifier 
# combinations with default parameters
################################################################################
library(randomForest)
library(plyr)
library(ROCR)
source("/hps/nobackup/ma/natalja/data/mQTL_data/scripts/load_data.R")

nr <- 500
res_set <- data.frame(nr = c(1:nr))
sets <-c("df1_a","df2_a","df3_a","df4_a",
"df1_b","df2_b","df3_b","df4_b",
"df1_c","df2_c","df3_c","df4_c",
"df1_d","df2_d","df3_d","df4_d")

for (set in sets){
	data_set <- eval(parse(text = set))
	err_rates <- c(0)
	for(i in 1:10){
		rf=randomForest(Diagnosis ~ ., data = data_set)
		err_rates <- rbind(err_rates,rf$err.rate[,1])
	}
	err_rates <- err_rates[-1,]
	res_set <- cbind(res_set,apply(err_rates,2,mean))
}

cols <- c("#f58231",
  "#f58231",
  "#f58231",
  "#f58231",
  "#000080",
  "#000080",
  "#000080",
  "#000080",
  "#808000",
  "#808000",
  "#808000",
  "#808000",
  "#e6194b",
  "#e6194b",
  "#e6194b",
  "#e6194b")

ltys <- c(2,3,4,1,2,3,4,1,2,3,4,1,2,3,4,1)

png("RF_all_sets.png")
matplot(1:nr , res_set, col=cols,type="l",lty=ltys, ylab="OOB error",xlab="trees")
legend("topright",
       legend=c("Set A, Classifier I","Set A, Classifier II","Set A, Classifier III","Set A, Classifier IV",
                "Set B, Classifier I","Set B, Classifier II","Set B, Classifier III","Set B, Classifier IV",
                "Set C, Classifier I","Set C, Classifier II","Set C, Classifier III","Set C, Classifier IV",
                "Set D, Classifier I","Set D, Classifier II","Set D, Classifier III","Set D, Classifier IV")
       ,lty=ltys, 
       col=cols,cex=0.8)
dev.off()


######################
nr <- 500
res_set <- data.frame(nr = c(1:nr))
sets <-c("df1_a","df2_a","df4_a",
         "df1_b","df2_b","df4_b",
         "df1_c","df2_c","df4_c",
         "df1_d","df2_d","df4_d")

for (set in sets){
  data_set <- eval(parse(text = set))
  err_rates <- c(0)
  for(i in 1:10){
    rf=randomForest(Diagnosis ~ ., data = data_set)
    err_rates <- rbind(err_rates,rf$err.rate[,1])
  }
  err_rates <- err_rates[-1,]
  res_set <- cbind(res_set,apply(err_rates,2,mean))
}
cols <- c("#f58231",
          "#000080",
          "#808000",
          "#f58231",
          "#000080",
          "#808000",
          "#f58231",
          "#000080",
          "#808000",
          "#f58231",
          "#000080",
          "#808000")

ltys <- c(4,4,4,3,3,3,2,2,2,1,1,1)

png("RF_all_sets3.png", width = 6, height = 6, units = 'in', res=300)
matplot(1:nr , res_set[,-1], col=cols,type="l",lty=ltys, ylab="OOB error",xlab="trees")
legend("topright",
       legend=c("Set A, Classifier I","Set A, Classifier II","Set A, Classifier III",
                "Set B, Classifier I","Set B, Classifier II","Set B, Classifier III",
                "Set C, Classifier I","Set C, Classifier II","Set C, Classifier III",
                "Set D, Classifier I","Set D, Classifier II","Set D, Classifier III")
       ,lty=ltys, 
       col=cols,cex=0.7)
dev.off()

#########################
# ROC curves
# D
sets <-c("df1_d","df2_d","df3_d","df4_d")
data <- df1_d
data_MDI <- data[data$Diagnosis %in% c("sMCI","cMCI"),]
for (set in sets){
  data <- eval(parse(text = set))
  data <- data[data$Diagnosis %in% c("ADC","CTL"),]
  data$ID <- rownames(data)
  data$Diagnosis <- relevel(data$Diagnosis,ref="CTL")
  data$Diagnosis <- as.factor(as.character(data$Diagnosis))
  test <- ddply(data,.(Diagnosis),function(x) x[sample(nrow(x),20),])
  train <- data[!data$ID %in% test$ID,]
  drops <- c("ID")
  train <- train [ , !(names(train ) %in% drops)]
  test <- test [ , !(names(test ) %in% drops)]
  
  # ADvsCTL values
  RF_train=randomForest(Diagnosis ~ ., data = train, importance = TRUE, proximity=TRUE)
  RF_predict = predict(RF_train, type = "prob", newdata = test)
  RF_prediction = prediction(RF_predict[,2], as.numeric(test$Diagnosis)-1)
  RF_performance = performance(RF_prediction, "tpr", "fpr")
  
  perf_AUC=performance(RF_prediction,"auc") #Calculate the AUC value
  AUC=perf_AUC@y.values[[1]]
  assign(paste('AUC_', parse(text = set),'_AD', sep=''), AUC) 
  assign(paste('RF_performance_', parse(text = set),'_AD', sep=''), RF_performance) 
  
  # cMCIvssMCI values
  RF_predict = predict(RF_train, type = "prob", newdata = data_MDI)
  RF_prediction = prediction(RF_predict[,2], as.numeric(data_MDI$Diagnosis)-1)
  table(data_MDI$Diagnosis, ifelse(RF_predict[,2] > 0.5,"CTL","ADC"))
  RF_performance = performance(RF_prediction, "tpr", "fpr")
  perf_AUC=performance(RF_prediction,"auc") #Calculate the AUC value
  AUC=perf_AUC@y.values[[1]]
  assign(paste('AUC_', parse(text = set),'_MDI', sep=''), AUC) 
  assign(paste('RF_performance_', parse(text = set),'_MDI', sep=''), RF_performance) 
}

# A
sets <-c("df1_a","df2_a","df3_a","df4_a")
data <- df1_a
data_MDI <- data[data$Diagnosis %in% c("sMCI","cMCI"),]
for (set in sets){
  data <- eval(parse(text = set))
  data <- data[data$Diagnosis %in% c("ADC","CTL"),]
  data$ID <- rownames(data)
  data$Diagnosis <- relevel(data$Diagnosis,ref="CTL")
  data$Diagnosis <- as.factor(as.character(data$Diagnosis))
  test <- ddply(data,.(Diagnosis),function(x) x[sample(nrow(x),20),])
  train <- data[!data$ID %in% test$ID,]
  drops <- c("ID")
  train <- train [ , !(names(train ) %in% drops)]
  test <- test [ , !(names(test ) %in% drops)]
  
  # ADvsCTL values
  RF_train=randomForest(Diagnosis ~ ., data = train, importance = TRUE, proximity=TRUE)
  RF_predict = predict(RF_train, type = "prob", newdata = test)
  RF_prediction = prediction(RF_predict[,2], as.numeric(test$Diagnosis)-1)
  RF_performance = performance(RF_prediction, "tpr", "fpr")
  perf_AUC=performance(RF_prediction,"auc") #Calculate the AUC value
  AUC=perf_AUC@y.values[[1]]
  assign(paste('AUC_', parse(text = set),'_AD', sep=''), AUC) 
  assign(paste('RF_performance_', parse(text = set),'_AD', sep=''), RF_performance) 
  
  # cMCIvssMCI values
  RF_predict = predict(RF_train, type = "prob", newdata = data_MDI)
  RF_prediction = prediction(RF_predict[,2], as.numeric(data_MDI$Diagnosis)-1)
  table(data_MDI$Diagnosis, ifelse(RF_predict[,2] > 0.5,"CTL","ADC"))
  RF_performance = performance(RF_prediction, "tpr", "fpr")
  perf_AUC=performance(RF_prediction,"auc") #Calculate the AUC value
  AUC=perf_AUC@y.values[[1]]
  assign(paste('AUC_', parse(text = set),'_MDI', sep=''), AUC) 
  assign(paste('RF_performance_', parse(text = set),'_MDI', sep=''), RF_performance) 
}

# B
sets <-c("df1_b","df2_b","df3_b","df4_b")
data <- df1_b
data_MDI <- data[data$Diagnosis %in% c("sMCI","cMCI"),]
for (set in sets){
  data <- eval(parse(text = set))
  data <- data[data$Diagnosis %in% c("ADC","CTL"),]
  data$ID <- rownames(data)
  data$Diagnosis <- relevel(data$Diagnosis,ref="CTL")
  data$Diagnosis <- as.factor(as.character(data$Diagnosis))
  test <- ddply(data,.(Diagnosis),function(x) x[sample(nrow(x),20),])
  train <- data[!data$ID %in% test$ID,]
  drops <- c("ID")
  train <- train [ , !(names(train ) %in% drops)]
  test <- test [ , !(names(test ) %in% drops)]
  
  # ADvsCTL values
  RF_train=randomForest(Diagnosis ~ ., data = train, importance = TRUE, proximity=TRUE)
  RF_predict = predict(RF_train, type = "prob", newdata = test)
  RF_prediction = prediction(RF_predict[,2], as.numeric(test$Diagnosis)-1)
  RF_performance = performance(RF_prediction, "tpr", "fpr")
  perf_AUC=performance(RF_prediction,"auc") #Calculate the AUC value
  AUC=perf_AUC@y.values[[1]]
  assign(paste('AUC_', parse(text = set),'_AD', sep=''), AUC) 
  assign(paste('RF_performance_', parse(text = set),'_AD', sep=''), RF_performance) 
  
  # cMCIvssMCI values
  RF_predict = predict(RF_train, type = "prob", newdata = data_MDI)
  RF_prediction = prediction(RF_predict[,2], as.numeric(data_MDI$Diagnosis)-1)
  table(data_MDI$Diagnosis, ifelse(RF_predict[,2] > 0.5,"CTL","ADC"))
  RF_performance = performance(RF_prediction, "tpr", "fpr")
  perf_AUC=performance(RF_prediction,"auc") #Calculate the AUC value
  AUC=perf_AUC@y.values[[1]]
  assign(paste('AUC_', parse(text = set),'_MDI', sep=''), AUC) 
  assign(paste('RF_performance_', parse(text = set),'_MDI', sep=''), RF_performance) 
}

# C
sets <-c("df1_c","df2_c","df3_c","df4_c")
data <- df1_c
data_MDI <- data[data$Diagnosis %in% c("sMCI","cMCI"),]
for (set in sets){
  data <- eval(parse(text = set))
  data <- data[data$Diagnosis %in% c("ADC","CTL"),]
  data$ID <- rownames(data)
  data$Diagnosis <- relevel(data$Diagnosis,ref="CTL")
  data$Diagnosis <- as.factor(as.character(data$Diagnosis))
  test <- ddply(data,.(Diagnosis),function(x) x[sample(nrow(x),20),])
  train <- data[!data$ID %in% test$ID,]
  drops <- c("ID")
  train <- train [ , !(names(train ) %in% drops)]
  test <- test [ , !(names(test ) %in% drops)]
  
  # ADvsCTL values
  RF_train=randomForest(Diagnosis ~ ., data = train, importance = TRUE, proximity=TRUE)
  RF_predict = predict(RF_train, type = "prob", newdata = test)
  RF_prediction = prediction(RF_predict[,2], as.numeric(test$Diagnosis)-1)
  RF_performance = performance(RF_prediction, "tpr", "fpr")
  
  perf_AUC=performance(RF_prediction,"auc") #Calculate the AUC value
  AUC=perf_AUC@y.values[[1]]
  assign(paste('AUC_', parse(text = set),'_AD', sep=''), AUC) 
  assign(paste('RF_performance_', parse(text = set),'_AD', sep=''), RF_performance) 
  
  # cMCIvssMCI values
  RF_predict = predict(RF_train, type = "prob", newdata = data_MDI)
  RF_prediction = prediction(RF_predict[,2], as.numeric(data_MDI$Diagnosis)-1)
  table(data_MDI$Diagnosis, ifelse(RF_predict[,2] > 0.5,"CTL","ADC"))
  RF_performance = performance(RF_prediction, "tpr", "fpr")
  perf_AUC=performance(RF_prediction,"auc") #Calculate the AUC value
  AUC=perf_AUC@y.values[[1]]
  assign(paste('AUC_', parse(text = set),'_MDI', sep=''), AUC) 
  assign(paste('RF_performance_', parse(text = set),'_MDI', sep=''), RF_performance) 
}


cols <- c("#f58231",
          "#f58231",
          "#f58231",
          "#f58231",
          "#000080",
          "#000080",
          "#000080",
          "#000080",
          "#808000",
          "#808000",
          "#808000",
          "#808000",
          "#e6194b",
          "#e6194b",
          "#e6194b",
          "#e6194b")

ltys <- c(2,3,4,1,2,3,4,1,2,3,4,1,2,3,4,1)

sets <-c("df1_a","df2_a","df3_a","df4_a",
         "df1_b","df2_b","df3_b","df4_b",
         "df1_c","df2_c","df3_c","df4_c",
         "df1_d","df2_d","df3_d","df4_d")
legends <- c("Set A, Classifier I","Set A, Classifier II","Set A, Classifier III","Set A, Classifier IV",
             "Set B, Classifier I","Set B, Classifier II","Set B, Classifier III","Set B, Classifier IV",
             "Set C, Classifier I","Set C, Classifier II","Set C, Classifier III","Set C, Classifier IV",
             "Set D, Classifier I","Set D, Classifier II","Set D, Classifier III","Set D, Classifier IV")


legends_AUC <- c()
i<-1
for (set in sets){
  AUC_val <- eval(parse(text=paste('AUC_', set,'_AD', sep='')))
  legends_AUC[i]=paste(legends[i],", AUC=",round(AUC_val,digits=3),sep="")
  i<-i+1
}

png("ROC_RF_A_sets_AD.png", width = 4, height = 4, units = 'in', res=300)
#par(mar=c(5.1,4.1,6,2.1), xpd=TRUE)
plot(RF_performance_df1_a_AD,col = cols[1],lty=ltys[1], lwd=2)
i<-2
for (set in sets[2:4]){
  prf <- eval(parse(text=paste('RF_performance_', set,'_AD', sep='')))
  plot(prf, add = TRUE, col = cols[i],lty=ltys[i], lwd=2)
  i <- i+1
}
abline(a=0, b= 1)
legend("bottomright",
       legend=legends_AUC[1:4],
       ,lty=ltys[1:4], 
       col=cols[1:4],cex=0.6)
dev.off()

png("ROC_RF_B_sets_AD.png", width = 4, height = 4, units = 'in', res=300)
#par(mar=c(5.1,4.1,6,2.1), xpd=TRUE)
plot(RF_performance_df1_b_AD,col = cols[5],lty=ltys[5], lwd=2)
i<-6
for (set in sets[6:8]){
  prf <- eval(parse(text=paste('RF_performance_', set,'_AD', sep='')))
  plot(prf, add = TRUE, col = cols[i],lty=ltys[i], lwd=2)
  i <- i+1
}
abline(a=0, b= 1)
legend("bottomright",
       legend=legends_AUC[5:8],
       ,lty=ltys[5:8], 
       col=cols[5:8],cex=0.6)
dev.off()


png("ROC_RF_C_sets_AD.png", width = 4, height = 4, units = 'in', res=300)
#par(mar=c(5.1,4.1,6,2.1), xpd=TRUE)
plot(RF_performance_df1_c_AD,col = cols[9],lty=ltys[9], lwd=2)
i<-10
for (set in sets[10:12]){
  prf <- eval(parse(text=paste('RF_performance_', set,'_AD', sep='')))
  plot(prf, add = TRUE, col = cols[i],lty=ltys[i], lwd=2)
  i <- i+1
}
abline(a=0, b= 1)
legend("bottomright",
       legend=legends_AUC[9:12],
       ,lty=ltys[9:12], 
       col=cols[9:12],cex=0.6)
dev.off()

png("ROC_RF_D_sets_AD.png", width = 4, height = 4, units = 'in', res=300)
#par(mar=c(5.1,4.1,6,2.1), xpd=TRUE)
plot(RF_performance_df1_d_AD,col = cols[13],lty=ltys[13], lwd=2)
i<-14
for (set in sets[14:16]){
  prf <- eval(parse(text=paste('RF_performance_', set,'_AD', sep='')))
  plot(prf, add = TRUE, col = cols[i],lty=ltys[i], lwd=2)
  i <- i+1
}
abline(a=0, b= 1)
legend("bottomright",
       legend=legends_AUC[13:16],
       ,lty=ltys[13:16], 
       col=cols[13:16],cex=0.6)
dev.off()

#MDI
png("ROC_RF_D_sets_MDI.png", width = 4, height = 4, units = 'in', res=300)
#par(mar=c(5.1,4.1,6,2.1), xpd=TRUE)
plot(RF_performance_df1_d_MDI,col = cols[13],lty=ltys[13], lwd=2)
i<-14
for (set in sets[14:16]){
  prf <- eval(parse(text=paste('RF_performance_', set,'_MDI', sep='')))
  plot(prf, add = TRUE, col = cols[i],lty=ltys[i], lwd=2)
  i <- i+1
}
abline(a=0, b= 1)
legend("bottomright",
       legend=legends_AUC[13:16],
       ,lty=ltys[13:16], 
       col=cols[13:16],cex=0.6)
dev.off()


