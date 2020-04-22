##################################################################################################
# Feature selection methods
# - Linear Regression with diagnosis in a model
# - WEKA
# - Random Forests
##################################################################################################

# Linear Regression with diagnosis in a model 
feature_selection_mqtl <- function(base){
  threshold <- 0.01
  dna_matrices_dir<-"./data/dna_matrices/"
  met_matrices_dir<-"./data/metabolomics_matrices/normalised/"
  cov_matrices_dir<-"./data/covariates/"
  res <- read.table(paste("./results/mqtl_",base,".txt",sep=""),header=T,sep="\t")
  res_sig <- res[res$FDR<threshold,]
  
  write.table(unique(res_sig$SNP),paste(base,"_SNP_sig.txt",sep=""),sep="\t",row.names=F,quote=F)
  write.table(unique(res_sig$gene),paste(base,"_met_sig.txt",sep=""),sep="\t",row.names=F,quote=F)
  
  # Matrices are to big to be processed in R
  # print(paste("grep -wFf ",base,"_sig_SNP.txt ",dna_matrices_dir,"dna_matrix_",base,".tsv > dna_",base,"_sig.txt",sep=""))
  # print(paste("grep -wFf ",base,"_sig_met.txt ",met_matrices_dir,base,".tsv > met_",base,"_sig.txt",sep=""))
  # print(paste("head -1 ",met_matrices_dir,base,".tsv > ",base,"_header.txt",sep=""))
  # print(paste("cat ",base,"_header.txt met_",base,"_sig.txt > met_",base,"_sig2.txt",sep=""))
  # print(paste("cat ",base,"_header.txt dna_",base,"_sig.txt > dna_",base,"_sig2.txt",sep=""))
   
}



################################################################################
# Random forests for feature selection
################################################################################
expr1 <- read.table("./results/selected_matrices/UHPOS_met_sig.txt",sep="\t",header=T,row.names=1,stringsAsFactors = F)
expr2 <- read.table("./results/selected_matrices/URPOS_met_sig.txt",sep="\t",header=T,row.names=1,stringsAsFactors = F)
expr3 <- read.table("./results/selected_matrices/URNEG_met_sig.txt",sep="\t",header=T,row.names=1,stringsAsFactors = F)
expr_nmr <- read.table("./results/selected_matrices/NMR_met_sig.txt",sep="\t",header=T,row.names=1,stringsAsFactors = F)

expr_nmr <- expr_nmr[,colnames(expr_nmr) %in% colnames(expr1)]
expr1 <- expr1[,colnames(expr1) %in% colnames(expr_nmr)]
expr2 <- expr2[,colnames(expr2) %in% colnames(expr_nmr)]
expr3 <- expr3[,colnames(expr3) %in% colnames(expr_nmr)]

expr_nmr <- expr_nmr[,colnames(expr1)]
expr2 <- expr2[,colnames(expr1)]
expr3 <- expr3[,colnames(expr1)]

expr <- rbind(expr_nmr,expr1)
expr <- rbind(expr,expr2)
expr <- rbind(expr,expr3)

cvrt = read.table("./data/covariates/Covariates_EigenMS_format.txt",sep="\t",header=T,row.names=1)
cvrt <- cvrt[rownames(cvrt) %in% colnames(expr),]
cvrt <- cvrt[colnames(expr),]

# a) Metabolites only and d) Metabolites and covariates
df <- as.data.frame(t(expr))
df$Gender <- cvrt$Gender
df$Centre <- cvrt$Centre
df$Age <- cvrt$Age
df$Diagnosis1 <- as.character(cvrt$Diagnosis)
df$Diagnosis2 <- as.character(cvrt$Diagnosis)
df$Diagnosis3 <- as.character(cvrt$Diagnosis)
df$Diagnosis2[df$Diagnosis2=="cMCI"] <- "ADC"
df$Diagnosis2[df$Diagnosis2=="sMCI"] <- "CTL"
df$Diagnosis3[df$Diagnosis3=="cMCI"] <- "ADC"
df$Diagnosis3[df$Diagnosis3=="sMCI"] <- "ADC"

# a) 
drops <- c("Gender","Centre","Age","Diagnosis2","Diagnosis3")
df1_a <- df[,!(names(df) %in% drops)]
names(df1_a)[names(df1_a) == 'Diagnosis1'] <- 'Diagnosis'

drops <- c("Gender","Centre","Age","Diagnosis1","Diagnosis3")
df2_a <- df[,!(names(df) %in% drops)]
names(df2_a)[names(df2_a) == 'Diagnosis2'] <- 'Diagnosis'

drops <- c("Gender","Centre","Age","Diagnosis1","Diagnosis2")
df3_a <- df[,!(names(df) %in% drops)]
names(df3_a)[names(df3_a) == 'Diagnosis3'] <- 'Diagnosis'

drops <- c("Gender","Centre","Age","Diagnosis2","Diagnosis3")
df4_a <- df[cvrt$Diagnosis %in% c("CTL","ADC"),!(names(df) %in% drops)]
names(df4_a)[names(df4_a) == 'Diagnosis1'] <- 'Diagnosis'

# or set d
drops <- c("Diagnosis2","Diagnosis3")
df1_d <- df[,!(names(df) %in% drops)]
names(df1_d)[names(df1_d) == 'Diagnosis1'] <- 'Diagnosis'

drops <- c("Diagnosis1","Diagnosis3")
df2_d <- df[,!(names(df) %in% drops)]
names(df2_d)[names(df2_d) == 'Diagnosis2'] <- 'Diagnosis'

drops <- c("Diagnosis1","Diagnosis2")
df3_d <- df[,!(names(df) %in% drops)]
names(df3_d)[names(df3_d) == 'Diagnosis3'] <- 'Diagnosis'

drops <- c("Diagnosis2","Diagnosis3")
df4_d <- df[cvrt$Diagnosis %in% c("CTL","ADC"),!(names(df) %in% drops)]
names(df4_d)[names(df4_d) == 'Diagnosis1'] <- 'Diagnosis'

# set b
dna_all <- read.table("./results/selected_matrices/dna_all.txt",sep="\t",header=T,row.names=1,stringsAsFactors = F)
# The common samples for NMR and LC-MS
dna_all<- dna_all[,colnames(dna_all) %in% colnames(expr)]
expr_dna <- expr[,colnames(expr) %in% colnames(dna_all)]
dna_all <- dna_all[,colnames(expr_dna)]
dna_all [is.na(dna_all )] <- 0
cvrt = read.table("./data/covariates/Covariates_EigenMS_format.txt",sep="\t",header=T,row.names=1)
cvrt <- cvrt[rownames(cvrt) %in% colnames(expr_dna),]
cvrt <- cvrt[colnames(expr_dna),]

df <- as.data.frame(cbind(t(expr_dna),t(dna_all)))
df$Gender <- cvrt$Gender
df$Centre <- cvrt$Centre
df$Age <- cvrt$Age
df$Diagnosis1 <- as.character(cvrt$Diagnosis)
df$Diagnosis2 <- as.character(cvrt$Diagnosis)
df$Diagnosis3 <- as.character(cvrt$Diagnosis)
df$Diagnosis2[df$Diagnosis2=="cMCI"] <- "ADC"
df$Diagnosis2[df$Diagnosis2=="sMCI"] <- "CTL"
df$Diagnosis3[df$Diagnosis3=="cMCI"] <- "ADC"
df$Diagnosis3[df$Diagnosis3=="sMCI"] <- "ADC"

# b) 
drops <- c("Gender","Centre","Age","Diagnosis2","Diagnosis3")
df1_b <- df[,!(names(df) %in% drops)]
names(df1_b)[names(df1_b) == 'Diagnosis1'] <- 'Diagnosis'

drops <- c("Gender","Centre","Age","Diagnosis1","Diagnosis3")
df2_b <- df[,!(names(df) %in% drops)]
names(df2_b)[names(df2_b) == 'Diagnosis2'] <- 'Diagnosis'

drops <- c("Gender","Centre","Age","Diagnosis1","Diagnosis2")
df3_b <- df[,!(names(df) %in% drops)]
names(df3_b)[names(df3_b) == 'Diagnosis3'] <- 'Diagnosis'

drops <- c("Gender","Centre","Age","Diagnosis2","Diagnosis3")
df4_b <- df[cvrt$Diagnosis %in% c("CTL","ADC"),!(names(df) %in% drops)]
names(df4_b)[names(df4_b) == 'Diagnosis1'] <- 'Diagnosis'

# or set c
drops <- c("Diagnosis2","Diagnosis3")
df1_c <- df[,!(names(df) %in% drops)]
names(df1_c)[names(df1_c) == 'Diagnosis1'] <- 'Diagnosis'

drops <- c("Diagnosis1","Diagnosis3")
df2_c <- df[,!(names(df) %in% drops)]
names(df2_c)[names(df2_c) == 'Diagnosis2'] <- 'Diagnosis'

drops <- c("Diagnosis1","Diagnosis2")
df3_c <- df[,!(names(df) %in% drops)]
names(df3_c)[names(df3_c) == 'Diagnosis3'] <- 'Diagnosis'

drops <- c("Diagnosis2","Diagnosis3")
df4_c <- df[cvrt$Diagnosis %in% c("CTL","ADC"),!(names(df) %in% drops)]
names(df4_c)[names(df4_c) == 'Diagnosis1'] <- 'Diagnosis'

df1_a$Diagnosis<- as.factor(df1_a$Diagnosis)
df1_b$Diagnosis<- as.factor(df1_b$Diagnosis)
df1_c$Diagnosis<- as.factor(df1_c$Diagnosis)
df1_d$Diagnosis<- as.factor(df1_d$Diagnosis)
df2_a$Diagnosis<- as.factor(df2_a$Diagnosis)
df2_b$Diagnosis<- as.factor(df2_b$Diagnosis)
df2_c$Diagnosis<- as.factor(df2_c$Diagnosis)
df2_d$Diagnosis<- as.factor(df2_d$Diagnosis)
df3_a$Diagnosis<- as.factor(df3_a$Diagnosis)
df3_b$Diagnosis<- as.factor(df3_b$Diagnosis)
df3_c$Diagnosis<- as.factor(df3_c$Diagnosis)
df3_d$Diagnosis<- as.factor(df3_d$Diagnosis)
df4_a$Diagnosis<- as.factor(df4_a$Diagnosis)
df4_b$Diagnosis<- as.factor(df4_b$Diagnosis)
df4_c$Diagnosis<- as.factor(df4_c$Diagnosis)
df4_d$Diagnosis<- as.factor(df4_d$Diagnosis)
names(df1_a) <- make.names(names(df1_a))
names(df1_b) <- make.names(names(df1_b))
names(df1_c) <- make.names(names(df1_c))
names(df1_d) <- make.names(names(df1_d))
names(df2_a) <- make.names(names(df2_a))
names(df2_b) <- make.names(names(df2_b))
names(df2_c) <- make.names(names(df2_c))
names(df2_d) <- make.names(names(df2_d))
names(df3_a) <- make.names(names(df3_a))
names(df3_b) <- make.names(names(df3_b))
names(df3_c) <- make.names(names(df3_c))
names(df3_d) <- make.names(names(df3_d))
names(df4_a) <- make.names(names(df4_a))
names(df4_b) <- make.names(names(df4_b))
names(df4_c) <- make.names(names(df4_c))
names(df4_d) <- make.names(names(df4_d)) 


nr <- 500
################################################################################
# Set A
############################
df_selected <- df1_a
names(df_selected) <- make.names(names(df_selected))
df_selected$Diagnosis <- as.factor(df_selected$Diagnosis)
samp_1_a <- sample(nrow(df_selected), 0.6 * nrow(df_selected))
test <- df_selected[samp_1_a,]
train <- df_selected[-samp_1_a,]
# dfeault number of trees 500, default mtry = 39
rf1_a=randomForest(Diagnosis ~ ., data = df_selected, ntree=nr)
df_selected <- df2_a
names(df_selected) <- make.names(names(df_selected))
df_selected$Diagnosis <- as.factor(df_selected$Diagnosis)
samp_2_a <- sample(nrow(df_selected), 0.6 * nrow(df_selected))
test <- df_selected[samp_2_a,]
train <- df_selected[-samp_2_a,]
# dfeault number of trees 500, default mtry = 39
rf2_a=randomForest(Diagnosis ~ ., data = df_selected, ntree=nr)
df_selected <- df3_a
names(df_selected) <- make.names(names(df_selected))
df_selected$Diagnosis <- as.factor(df_selected$Diagnosis)
samp_3_a <- sample(nrow(df_selected), 0.6 * nrow(df_selected))
test <- df_selected[samp_3_a,]
train <- df_selected[-samp_3_a,]
# dfeault number of trees 500, default mtry = 39
rf3_a=randomForest(Diagnosis ~ ., data = df_selected, ntree=nr)
df_selected <- df4_a
names(df_selected) <- make.names(names(df_selected))
df_selected$Diagnosis <- as.factor(df_selected$Diagnosis)
samp_4_a <- sample(nrow(df_selected), 0.6 * nrow(df_selected))
test <- df_selected[samp_4_a,]
train <- df_selected[-samp_4_a,]
# dfeault number of trees 500, default mtry = 39
rf4_a=randomForest(Diagnosis ~ ., data = df_selected, ntree=nr)

################################################################################
# Set B
############################
df_selected <- df1_b
names(df_selected) <- make.names(names(df_selected))
df_selected$Diagnosis <- as.factor(df_selected$Diagnosis)
samp_1_b <- sample(nrow(df_selected), 0.6 * nrow(df_selected))
test <- df_selected[samp_1_b,]
train <- df_selected[-samp_1_b,]
# dfeault number of trees 500, default mtry = 39
rf1_b=randomForest(Diagnosis ~ ., data = df_selected, ntree=nr)
df_selected <- df2_b
names(df_selected) <- make.names(names(df_selected))
df_selected$Diagnosis <- as.factor(df_selected$Diagnosis)
samp_2_b <- sample(nrow(df_selected), 0.6 * nrow(df_selected))
test <- df_selected[samp_2_b,]
train <- df_selected[-samp_2_b,]
# dfeault number of trees 500, default mtry = 39
rf2_b=randomForest(Diagnosis ~ ., data = df_selected, ntree=nr)
df_selected <- df3_b
names(df_selected) <- make.names(names(df_selected))
df_selected$Diagnosis <- as.factor(df_selected$Diagnosis)
samp_3_b <- sample(nrow(df_selected), 0.6 * nrow(df_selected))
test <- df_selected[samp_3_b,]
train <- df_selected[-samp_3_b,]
# dfeault number of trees 500, default mtry = 39
rf3_b=randomForest(Diagnosis ~ ., data = df_selected, ntree=nr)
df_selected <- df4_b
names(df_selected) <- make.names(names(df_selected))
df_selected$Diagnosis <- as.factor(df_selected$Diagnosis)
samp_4_b <- sample(nrow(df_selected), 0.6 * nrow(df_selected))
test <- df_selected[samp_4_b,]
train <- df_selected[-samp_4_b,]
# dfeault number of trees 500, default mtry = 39
rf4_b=randomForest(Diagnosis ~ ., data = df_selected, ntree=nr)

################################################################################
# Set C
############################
df_selected <- df1_c
names(df_selected) <- make.names(names(df_selected))
df_selected$Diagnosis <- as.factor(df_selected$Diagnosis)
samp_1_c <- sample(nrow(df_selected), 0.6 * nrow(df_selected))
test <- df_selected[samp_1_c,]
train <- df_selected[-samp_1_c,]
# dfeault number of trees 500, default mtry = 39
rf1_c=randomForest(Diagnosis ~ ., data = df_selected, ntree=nr)
df_selected <- df2_c
names(df_selected) <- make.names(names(df_selected))
df_selected$Diagnosis <- as.factor(df_selected$Diagnosis)
samp_2_c <- sample(nrow(df_selected), 0.6 * nrow(df_selected))
test <- df_selected[samp_2_c,]
train <- df_selected[-samp_2_c,]
# dfeault number of trees 500, default mtry = 39
rf2_c=randomForest(Diagnosis ~ ., data = df_selected, ntree=nr)
df_selected <- df3_c
names(df_selected) <- make.names(names(df_selected))
df_selected$Diagnosis <- as.factor(df_selected$Diagnosis)
samp_3_c <- sample(nrow(df_selected), 0.6 * nrow(df_selected))
test <- df_selected[samp_3_c,]
train <- df_selected[-samp_3_c,]
# dfeault number of trees 500, default mtry = 39
rf3_c=randomForest(Diagnosis ~ ., data = df_selected, ntree=nr)
df_selected <- df4_c
names(df_selected) <- make.names(names(df_selected))
df_selected$Diagnosis <- as.factor(df_selected$Diagnosis)
samp_4_c <- sample(nrow(df_selected), 0.6 * nrow(df_selected))
test <- df_selected[samp_4_c,]
train <- df_selected[-samp_4_c,]
# dfeault number of trees 500, default mtry = 39
rf4_c=randomForest(Diagnosis ~ ., data = df_selected, ntree=nr)

################################################################################
# Set D
############################
df_selected <- df1_d
names(df_selected) <- make.names(names(df_selected))
df_selected$Diagnosis <- as.factor(df_selected$Diagnosis)
samp_1_d <- sample(nrow(df_selected), 0.6 * nrow(df_selected))
test <- df_selected[samp_1_d,]
train <- df_selected[-samp_1_d,]
# dfeault number of trees 500, default mtry = 39
rf1_d=randomForest(Diagnosis ~ ., data = df_selected, ntree=nr)
df_selected <- df2_d
names(df_selected) <- make.names(names(df_selected))
df_selected$Diagnosis <- as.factor(df_selected$Diagnosis)
samp_2_d <- sample(nrow(df_selected), 0.6 * nrow(df_selected))
test <- df_selected[samp_2_d,]
train <- df_selected[-samp_2_d,]
# dfeault number of trees 500, default mtry = 39
rf2_d=randomForest(Diagnosis ~ ., data = df_selected, ntree=nr)
df_selected <- df3_d
names(df_selected) <- make.names(names(df_selected))
df_selected$Diagnosis <- as.factor(df_selected$Diagnosis)
samp_3_d <- sample(nrow(df_selected), 0.6 * nrow(df_selected))
test <- df_selected[samp_3_d,]
train <- df_selected[-samp_3_d,]
# dfeault number of trees 500, default mtry = 39
rf3_d=randomForest(Diagnosis ~ ., data = df_selected, ntree=nr)
df_selected <- df4_d
names(df_selected) <- make.names(names(df_selected))
df_selected$Diagnosis <- as.factor(df_selected$Diagnosis)
samp_4_d <- sample(nrow(df_selected), 0.6 * nrow(df_selected))
test <- df_selected[samp_4_d,]
train <- df_selected[-samp_4_d,]
# dfeault number of trees 500, default mtry = 39
rf4_d=randomForest(Diagnosis ~ ., data = df_selected, ntree=nr)

library(randomForest)
source("/hps/nobackup/ma/natalja/data/mQTL_data/scripts/load_data.R")

err_rates <- c(0)
for(i in 1:10){
  rf=randomForest(Diagnosis ~ ., data = df4_b)
  err_rates <- rbind(err_rates,rf$err.rate[,1])
}
err_rates <- err_rates[-1,]
mean_res4_b <- apply(err_rates,2,mean)

err_rates <- c(0)
for(i in 1:10){
  rf=randomForest(Diagnosis ~ ., data = df4_a)
  err_rates <- rbind(err_rates,rf$err.rate[,1])
}
err_rates <- err_rates[-1,]
mean_res4_a <- apply(err_rates,2,mean)

err_rates <- c(0)
for(i in 1:10){
  rf=randomForest(Diagnosis ~ ., data = df4_c)
  err_rates <- rbind(err_rates,rf$err.rate[,1])
}
err_rates <- err_rates[-1,]
mean_res4_c <- apply(err_rates,2,mean)

err_rates <- c(0)
for(i in 1:10){
  rf=randomForest(Diagnosis ~ ., data = df4_d)
  err_rates <- rbind(err_rates,rf$err.rate[,1])
}
err_rates <- err_rates[-1,]
mean_res4_d <- apply(err_rates,2,mean)


mean_B4 <- mean_res
mean_ADC <- mean_res
err_rates <- cbind(mean_ADC,mean_B4)

cols <- c("#f58231",
          "#000080",
          "#e6194b"
)

ltys <- c(1,1,1)
nr<-500
pdf("RF_ADC_B4_B.pdf")
matplot(1:nr , err_rates, col=cols,type="l",lty=ltys, ylab="Error",xlab="trees")
legend("topright",
       legend=c("Set A with 9 GWAS SNPs, Classifier IV","Set A with B4 SNPs, Classifier IV","Set B, Classifier IV")
       ,lty=ltys, 
       col=cols,cex=0.5)
dev.off()


rf1_a <- randomForest(Diagnosis ~ ., data = df1_a)
rf2_a <- randomForest(Diagnosis ~ ., data = df2_a)
rf3_a <- randomForest(Diagnosis ~ ., data = df3_a)

rf1_b <- randomForest(Diagnosis ~ ., data = df1_b)
rf2_b <- randomForest(Diagnosis ~ ., data = df2_b)
rf3_b <- randomForest(Diagnosis ~ ., data = df3_b)

rf1_c <- randomForest(Diagnosis ~ ., data = df1_c)
rf2_c <- randomForest(Diagnosis ~ ., data = df2_c)
rf3_c <- randomForest(Diagnosis ~ ., data = df3_c)

rf1_d <- randomForest(Diagnosis ~ ., data = df1_d)
rf2_d <- randomForest(Diagnosis ~ ., data = df2_d)
rf3_d <- randomForest(Diagnosis ~ ., data = df3_d)

rf4_a <- randomForest(Diagnosis ~ ., data = df4_a)
rf4_b <- randomForest(Diagnosis ~ ., data = df4_b)
rf4_c <- randomForest(Diagnosis ~ ., data = df4_c)
rf4_d <- randomForest(Diagnosis ~ ., data = df4_d)

err_rates <- cbind(rf1_a$err.rate[,1],rf2_a$err.rate[,1])
err_rates <- cbind(err_rates,rf3_a$err.rate[,1])
err_rates <- cbind(err_rates,mean_res4_a)
err_rates <- cbind(err_rates,rf1_b$err.rate[,1])
err_rates <- cbind(err_rates,rf2_b$err.rate[,1])
err_rates <- cbind(err_rates,rf3_b$err.rate[,1])
err_rates <- cbind(err_rates,mean_res4_b)
err_rates <- cbind(err_rates,rf1_c$err.rate[,1])
err_rates <- cbind(err_rates,rf2_c$err.rate[,1])
err_rates <- cbind(err_rates,rf3_c$err.rate[,1])
err_rates <- cbind(err_rates,mean_res4_c)
err_rates <- cbind(err_rates,rf1_d$err.rate[,1])
err_rates <- cbind(err_rates,rf2_d$err.rate[,1])
err_rates <- cbind(err_rates,rf3_d$err.rate[,1])
err_rates <- cbind(err_rates,mean_res4_d)

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
matplot(1:nr , err_rates, col=cols,type="l",lty=ltys, ylab="OOB error",xlab="trees")
legend("topright",
       legend=c("Set A, Classifier I","Set A, Classifier II","Set A, Classifier III","Set A, Classifier IV",
                "Set B, Classifier I","Set B, Classifier II","Set B, Classifier III","Set B, Classifier IV",
                "Set C, Classifier I","Set C, Classifier II","Set C, Classifier III","Set C, Classifier IV",
                "Set D, Classifier I","Set D, Classifier II","Set D, Classifier III","Set D, Classifier IV")
       ,lty=ltys, 
       col=cols,cex=0.8)
dev.off()

min(err_rates)
#...
# Set B Diagnosis IV - smallest error = 0.03125
which(rf4_b$err.rate[,1]==min(rf4_b$err.rate[,1]))
#309



df_selected <- df4_a
names(df_selected) <- make.names(names(df_selected))
df_selected$Diagnosis <- as.factor(df_selected$Diagnosis)
test <- df_selected[samp_4_b,]
train <- df_selected[-samp_4_b,]

# default number of trees 500, default mtry = 39
rf=randomForest(Diagnosis ~ ., data = train)

which(Diagnosis.rf$err.rate[,1]==min(Diagnosis.rf$err.rate[,1]))
# 116 149 150

#Plotting the Error vs Number of Trees Graph.
pdf("RF_Errors_vs_NumberofTrees_B4.pdf")
plot(rf4_b, main="Set B, diagnosis IV")
legend("topright",
       legend=c("ADC class",
                "CTL class",
                "Out of Bag")
       ,lty=c(2,3,1), 
       col=c("red","green","black"))
dev.off()

nr <-400

oob.err=double(nr)
test.err_in=double(nr)
test.err_out=double(nr)

#mtry is no of Variables randomly chosen at each split
for(mtry in  1:nr) 
{
    rf=randomForest(Diagnosis ~ . , data = df_selected,mtry=mtry,ntree=309) 
    oob.err[mtry] = rf$err.rate[309,1] #Error of all Trees fitted
   
    pred<-predict(rf,train) #Predictions on Test Set for each Tree
    res <- table(pred,train$Diagnosis)
    test.err_in[mtry]= 1 - (res[1,1]+res[2,2])/sum(res)

  
    pred<-predict(rf,test) #Predictions on Test Set for each Tree
    res <- table(pred,test$Diagnosis)
    test.err_out[mtry]= 1 - (res[1,1]+res[2,2])/sum(res)
  
    cat(mtry," ") #printing the output to the console
}

res <- cbind(oob.err,test.err_in,test.err_out)
res <- res[which(rowSums(res) > 0),] 



which(oob.err==min(oob.err))
#Plotting the Error vs Number of Trees Graph.
pdf("RF_OutofBagError_B4.pdf")
matplot(1:nr , oob.err, type="l",lty=1, col=c("red"),ylab="Mean Squared Error",xlab="Number of Predictors Considered at each Split")
#legend("topright",legend=c("Out of Bag Error"),pch=19, col=c("red"))
dev.off()
#A4 mtry = 91, 93
#B4 mtry = 199, 322

# rf=randomForest(Diagnosis ~ . , data = df_selected,ntree=442, mtry=91)
#randomForest(formula = Diagnosis ~ ., data = df_selected, ntree = 500,      mtry = 91) 
#Type of random forest: classification
#Number of trees: 500
#No. of variables tried at each split: 91

#OOB estimate of  error rate: 2.37%
#Confusion matrix:
#  ADC CTL class.error
#ADC 160   2  0.01234568
#CTL   5 128  0.03759398



#Number of trees: 150 (ntree)
#No. of variables tried at each split: 159 (mtry)


rf=randomForest(Diagnosis ~ . , data = train ,mtry=159,ntree=150,importance=TRUE) 
# oob.err = 2.82%
pred <- predict(rf, newdata = test)

table(pred, test$Diagnosis)

importance    <- importance(rf,type=1,scale=FALSE)
varImportance <- data.frame(Variables = row.names(importance), MDA = round(importance[ ,'MeanDecreaseAccuracy'],5), MDG = round(importance[ ,'MeanDecreaseGini'],5))

write.csv(varImportance,"feature_selected_D4.txt", row.names = FALSE)

importance_cfa_conditional <- varimp(cfa,conditional=TRUE)

varImportance_cfa <- data.frame(Variables = names(importance_cfa), CF = round(as.numeric(importance_cfa),5), 
                                CF_COND=round(as.numeric(importance_cfa_conditional),5))
write.csv(varImportance_cfa,"feature_selected_D4_cf.txt", row.names = FALSE)

#Importance scores
df_selected <- df4_b
names(df_selected) <- make.names(names(df_selected))
df_selected$Diagnosis <- as.factor(df_selected$Diagnosis)
df_selected_b <- df_selected

df_selected <- df4_a
names(df_selected) <- make.names(names(df_selected))
df_selected$Diagnosis <- as.factor(df_selected$Diagnosis)
df_selected_a <- df_selected

# mtry and ntree - parameters tuning results
rf_a=randomForest(Diagnosis ~ . , data = df_selected_a ,mtry=91,ntree=443,importance=TRUE)
rf_b=randomForest(Diagnosis ~ . , data = df_selected_b ,mtry=199,ntree=309,importance=TRUE)

importance    <- importance(rf_b)
varImportance <- data.frame(Variables = row.names(importance), MDA_b = round(importance[ ,'MeanDecreaseAccuracy'],5), MDG_b = round(importance[ ,'MeanDecreaseGini'],5))
varImportance$MDA_a <- 0
varImportance$MDG_a <- 0
importance    <- importance(rf_a)
varImportance_a <- data.frame(Variables = row.names(importance), MDA_a = round(importance[ ,'MeanDecreaseAccuracy'],5), MDG_a = round(importance[ ,'MeanDecreaseGini'],5))

for(v in varImportance_a$Variables){
  varImportance[varImportance$Variables==v,]$MDA_a <- varImportance_a[varImportance_a$Variables==v,]$MDA_a
  varImportance[varImportance$Variables==v,]$MDG_a <- varImportance_a[varImportance_a$Variables==v,]$MDG_a
}

library(party)
cfa <- cforest(Diagnosis ~ . ,data=df_selected_a,control=cforest_unbiased(mtry=91,ntree=443))

library(party)
cfa <- cforest(Diagnosis ~ . ,data=df4_d,control=cforest_unbiased(mtry=70,ntree=680))

importance_cfa <- varimp(cfa)

varImportance$CF_a <- 0
varImportance$CFCOND_a <- 0
varImportance$CF_b <- 0
varImportance$CFCOND_b <- 0
importance_cfa_conditional <- varimp(cfa,conditional=TRUE)

varImportance_cfa <- data.frame(Variables = names(importance_cfa), CF = round(as.numeric(importance_cfa),5), 
                                CF_COND=round(as.numeric(importance_cfa_conditional),5))

for(v in names(importance_cfa)){
  varImportance[varImportance$Variables==v,]$CF_a <-importance_cfa[names(importance_cfa)==v,]
  varImportance[varImportance$Variables==v,]$MDG_a <- varImportance_a[varImportance_a$Variables==v,]$MDG_a
}

cfb <- cforest(Diagnosis ~ . ,data=df_selected_b,control=cforest_unbiased(mtry=199,ntree=309))
importance_cfb <- varimp(cfb)
importance_cfb_conditional <- varimp(cfb,conditional=TRUE)
varImportance_cfb <- data.frame(Variables = names(importance_cfb), CF = round(as.numeric(importance_cfb),5), 
                                CF_COND=round(as.numeric(importance_cfb_conditional),5))

for(v in names(importance_cfb)){
  varImportance[varImportance$Variables==v,]$CF_a <-varImportance_cfa[varImportance_cfa$Variables==v,]$CF
  varImportance[varImportance$Variables==v,]$CFCOND_a <- varImportance_cfa[varImportance_cfa$Variables==v,]$CF_COND
  varImportance[varImportance$Variables==v,]$CF_b <-varImportance_cfb[varImportance_cfb$Variables==v,]$CF
  varImportance[varImportance$Variables==v,]$CFCOND_b <- varImportance_cfb[varImportance_cfb$Variables==v,]$CF_COND
}
# Recursive Feature Elimination incorporating resampling
# https://topepo.github.io/caret/recursive-feature-elimination.html#rfe
library(caret)
# define the control using a random forest selection function for 10 fold cross-validation
control<-rfeControl(functions=rfFuncs, method="cv", number=10)
dim(df_selected)
# run the RFE algorithm
results <- rfe(df_selected[,1:1542],df_selected[,1543], sizes=c(1:1542),rfeControl=control)
# summarize the results
print(results)
# list the chosen features
predictors(results)
# plot the results
pdf("RF_RFE_A4.pdf")
plot(results, type=c("g", "o"))
dev.off()
write.csv(results,"feature_selected_RFE_A4.txt", row.names = FALSE)


df_selected <- df4_b
names(df_selected) <- make.names(names(df_selected))
df_selected$Diagnosis <- as.factor(df_selected$Diagnosis)
rf=randomForest(Diagnosis ~ . , data = df_selected) 
df_selected <- df1_a
names(df_selected) <- make.names(names(df_selected))
df_selected$Diagnosis <- as.factor(df_selected$Diagnosis)

test <- df_selected[df_selected$Diagnosis %in% c("sMCI","cMCI"),]
test$pred <- predict(rf,test)

table(test[,8475], test[,8476])

# B4 model, OOB 2.94% 
#      ADC CTL
#cMCI  11  12
#sMCI  12  68



# 48% of cMCI samples are classified as ADC
# 85% of sMCI samples are classified as CTL

df_selected <- df4_a
names(df_selected) <- make.names(names(df_selected))
df_selected$Diagnosis <- as.factor(df_selected$Diagnosis)
rf=randomForest(Diagnosis ~ . , data = df_selected) 

# A4 model, OOB 2.71%
#      ADC CTL
#cMCI  16   7
#sMCI  17  63

#     ADC CTL
#cMCI  26  10
#sMCI  34 106


# 72% of cMCI samples are classified as ADC
# 75% of sMCI samples are classified as CTL

library(Boruta)
df_selected <- df4_a
# df_selected <- df4_b
names(df_selected) <- make.names(names(df_selected))
df_selected$Diagnosis <- as.factor(df_selected$Diagnosis)
# for B4 - genotypes 
# to treat them as factor data type
convert <- c(1543:8474)
df_selected[,convert] <- data.frame(apply(df_selected[convert], 2, as.factor)) 


boruta.train <- Boruta(Diagnosis~., data = df_selected, doTrace = 2)

pdf("Boruta_B4.pdf")
plot(boruta.train, xlab = "", xaxt = "n")
dev.off()
# The tentative attributes will be classified as confirmed or rejected 
# by comparing the median Z score of the attributes with the median Z score of the best shadow attribute
final.boruta <- TentativeRoughFix(boruta.train)
print(final.boruta)
getSelectedAttributes(final.boruta, withTentative = F)

varImportance$Boruta1 <- 0
varImportance$Boruta2 <- 0
varImportance[varImportance$Variables %in% getSelectedAttributes(boruta.train, withTentative = F),]$Boruta1 <- 1
varImportance[varImportance$Variables %in% getSelectedAttributes(boruta.train, withTentative = T),]$Boruta2 <- 1

boruta.df <- attStats(final.boruta)
print(boruta.df)

write.csv(getSelectedAttributes(boruta.train, withTentative = T),"Boruta_feature_selected_A4.txt", row.names = FALSE)


library(mlr) # based on weka?
df_selected <- df4_a
names(df_selected) <- make.names(names(df_selected))
df_selected$Diagnosis <- as.factor(df_selected$Diagnosis)
df_task <- makeClassifTask(data = df_selected, target = "Diagnosis")
fv2 = generateFilterValuesData(df_task, method = c("information.gain", "chi.squared"))


pdf("FilterMethods_A4.pdf")
plotFilterValues(fv2)
dev.off()

library(caret)
# In order to remove redundant features the caret R package has been used. 
# It allows to analyze a correlation matrix of features (metabolites, genotypes 
# and covariates ???) and found the attributes that can be removed.  
# correlation of features except diagnoses 
drops <- c("Diagnosis","Diagnosis3","Diagnosis4")
df_selected <- df[,!(names(df) %in% drops)]
df_selected$Gender <- as.numeric(df_selected$Gender)
df_selected$Centre <- as.numeric(df_selected$Centre)

correlationMatrix <- cor(df_selected)
highlyCorrelated <- findCorrelation(correlationMatrix, cutoff=0.75)
filteredMatrix <- df_selected[,-highlyCorrelated]

df_selected <- cbind(filteredMatrix,df[,8475])
colnames(df_selected)[ncol(df_selected)] <- "Diagnosis"

names(df_selected) <- make.names(names(df_selected))
df_selected$Diagnosis <- as.factor(df_selected$Diagnosis)
df_selected$Diagnosis2 <- as.factor(df_selected$Diagnosis2)





control <- trainControl(method="repeatedcv", number=10, repeats=3)
# train the model
model <- train(Diagnosis ~ . , data=df_selected, method="lvq", preProcess="scale", trControl=control)
# estimate variable importance
importance <- varImp(model, scale=FALSE)
# summarize importance
print(importance)
# plot importance
plot(importance)

control <- rfeControl(functions=rfFuncs, method="cv", number=10)
# Recursive Feature Elimination algorithm
# A Random Forest algorithm is used on each iteration to evaluate the model. 
# The algorithm is configured to explore all possible subsets of the attributes. 
results <- rfe(df_selected[,1:ncol(df_selected)], df_selected[,ncol(df_selected)], sizes=c(1:ncol(df_selected)), rfeControl=control)


library(randomForest)

df_selected <- df[,colnames(df) %in% selected_mets$V1 | colnames(df) %in% c("Diagnosis2","Diagnosis") ]

df_selected <- df

df_selected <- df[cvrt$Diagnosis %in% c("CTL","ADC"),]

names(df_selected) <- make.names(names(df_selected))
df_selected$Diagnosis <- as.factor(df_selected$Diagnosis)
df_selected$Diagnosis2 <- as.factor(df_selected$Diagnosis2)
df_selected$Gender <- as.factor(df_selected$Gender)
df_selected$Centre <- as.factor(df_selected$Centre)

samp <- sample(nrow(df_selected), 0.6 * nrow(df_selected))
test <- df_selected[samp,]
train <- df_selected[-samp,]



drops <- c("Gender","Centre")
df_selected <- df_selected[,!(names(df_selected) %in% drops)]

samp <- sample(nrow(df_selected), 0.6 * nrow(df_selected))
test <- df_selected[samp,]
train <- df_selected[-samp,]

model <- randomForest(Diagnosis2 ~ . - Diagnosis, data = train, importance=TRUE,ntree=500)
pred <- predict(model, newdata = test)

table(pred, test$Diagnosis2)
png("RandomForest.png")
varImpPlot(model)
dev.off()

model <- randomForest(Diagnosis2 ~ . - Diagnosis, data = train, importance=TRUE,ntree=2000)
importance    <- importance(model)
varImportance <- data.frame(Variables = row.names(importance), Importance = round(importance[ ,'MeanDecreaseGini'],2))

#Create a rank variable based on importance
library(dplyr)
rankImportance <- varImportance %>% mutate(Rank = paste0('#',dense_rank(desc(Importance))))

#Use ggplot2 to visualize the relative importance of variables
library(ggplot2)
library(ggthemes)
ggplot(rankImportance, aes(x = reorder(Variables, Importance), 
                           y = Importance, fill = Importance)) +
  geom_bar(stat='identity') + 
  geom_text(aes(x = Variables, y = 0.5, label = Rank),
            hjust=0, vjust=0.55, size = 4, colour = 'red') +
  labs(x = 'Variables') +
  coord_flip() + 
  theme_few()

library(party)
fit <- cforest(Diagnosis ~ . - Diagnosis2,data = train,controls=cforest_unbiased(ntree=2000, mtry=3))
fit2 <- cforest(Diagnosis2 ~ . - Diagnosis,data = train,controls=cforest_unbiased(ntree=2000, mtry=3))
pred <-  predict(fit, test, OOB=TRUE, type = "response")
table(pred,test$Diagnosis)
pred2 <-  predict(fit2, test, OOB=TRUE, type = "response")
table(pred2,test$Diagnosis2)

png("RandomForest2.png")
varImpPlot(fit)
dev.off()


#########
#clusters
#A fundamental question is how to determine the value of the parameter k. 
#If we looks at the percentage of variance explained as a function of the number of clusters: 
#One should choose a number of clusters so that adding another cluster doesn’t give much better modeling of the data. 
#More precisely, if one plots the percentage of variance explained by the clusters against the number of clusters, 
#the first clusters will add much information (explain a lot of variance), but at some point the marginal gain will drop, 
# giving an angle in the graph. The number of clusters is chosen at this point, hence the “elbow criterion”.
wssplot <- function(data, nc=500, seed=1234){
  wss <- (nrow(data)-1)*sum(apply(data,2,var))
  for (i in 2:nc){
    set.seed(seed)
    wss[i] <- sum(kmeans(data, centers=i)$withinss)}
  plot(1:nc, wss, type="b", xlab="Number of Clusters",
       ylab="Within groups sum of squares")}

png("WSSplot_NMR.png")
wssplot(expr_nmr)
dev.off()

expr_lcms <- rbind(expr1,expr2)
expr_lcms <- rbind(expr_lcms,expr3)

png("WSSplot_LCMS.png")
wssplot(expr_lcms)
dev.off()

cvrt = read.table("/hps/nobackup/ma/natalja/data/mQTL_data/covariates/Covariates_all_vertical_text2.txt",sep="\t",header=T,row.names=1)
cvrt <- cvrt[rownames(cvrt) %in% colnames(expr_nmr),]
cvrt <- cvrt[colnames(expr_nmr),]


cor.mtest <- function(mat, ...) {
  mat <- as.matrix(mat)
  n <- ncol(mat)
  p.mat<- matrix(NA, n, n)
  diag(p.mat) <- 0
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      tmp <- cor.test(mat[, i], mat[, j], ...)
      p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
    }
  }
  colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
  p.mat
}


mat_data <- data.matrix(expr_nmr)
rownames(mat_data) <- gsub("X\\.","",rownames(mat_data))
rownames(mat_data) <- gsub("X","",rownames(mat_data))
M <- cor(t(mat_data))
p.mat <- cor.mtest(t(mat_data))
# matrix of the p-value of the correlation

head(p.mat[, 1:5])



png("Corr_plot2.png")
corrplot(M, type="upper", order="hclust", 
         p.mat = p.mat, sig.level = 0.01, insig = "blank",addgrid.col = NA)

dev.off()


################################################################################
#First set the mtry to the default value (sqrt of total number of all predictors) 
#and search for the optimal ntree value. 
#To find the number of trees that correspond to a stable classifier, 
#we build random forest with different ntree values (100, 200, 300….,1,000). 
#We build 10 RF classifiers for each ntree value, record the OOB error rate and 
#see the number of trees where the out of bag error rate stabilizes and reach minimum.
library(caret)
library(party)
source("/hps/nobackup/ma/natalja/data/mQTL_data/scripts/load_data.R")
tuning_results <- data.frame(set="",nr=0, ntree=0, kappa=0,accuracy=0)
control <- trainControl(method="oob")
ntree_seq <- c(seq(0, 2500, by = 100)[-1])
sets <-c("df4_b") #c("df1_d","df2_d","df3_d","df4_d")
# c("df1_b","df2_b","df3_b","df4_b")
# "df1_c","df2_c","df3_c","df4_c",
# "df1_d","df2_d","df3_d","df4_d")
#"df1_a","df2_a","df3_a","df4_a",
for (set in sets){
  tuning_results <- data.frame(set="",nr=0, ntree=0, kappa=0,accuracy=0)
  for (ntree in ntree_seq) {
    data_set <- eval(parse(text = set))
    mtry_fixed <- sqrt(ncol(data_set))
    for(i in 1:10){
      mod <- train(Diagnosis ~ .,
                   data = data_set,
                   method = "cforest",
                   tuneGrid = data.frame(.mtry = mtry_fixed),
                   #trControl = control,
                   controls = cforest_unbiased(ntree = ntree))
      mod_res <- data.frame(set=set,nr=i,ntree=ntree,kappa=mod$results$Kappa, accuracy=mod$results$Accuracy)
      tuning_results <- rbind(tuning_results,mod_res)
    }

    print(ntree)
  }
  print(set)
  write.csv(tuning_results,paste(set,"_results.csv",sep=""),row.names = F)
} 

library(randomForest)
source("./scripts/load_data.R")
tuning_results <- data.frame(set="",nr=0, ntree=0, oob=0)
control <- trainControl(method="oob")
ntree_seq <- c(seq(0, 1000, by = 100)[-1])
sets <- #c("df1_d","df2_d","df3_d","df4_d")
 c("df4_d")
for (set in sets){
  tuning_results <- data.frame(set="",nr=0, ntree=0, oob=0)
  for (ntree in ntree_seq) {
    data_set <- eval(parse(text = set))
    mtry_fixed <- sqrt(ncol(data_set))
    for(i in 1:10){
      rf <- randomForest(Diagnosis ~ ., data = data_set, ntree=ntree)
      mod_res <- data.frame(set=set,nr=i,ntree=ntree,oob=rf$err.rate[nrow(rf$err.rate),1])
      tuning_results <- rbind(tuning_results,mod_res)
    }
    
    print(ntree)
  }
  print(set)
  write.csv(tuning_results,paste(set,"_results_rf.csv",sep=""),row.names = F)
} 


set <- "df4_d"
ntree_seq <- c(seq(200, 399, by = 10)[-1])
tuning_results <- data.frame(set="",nr=0, ntree=0, oob=0)
for (ntree in ntree_seq) {
  data_set <- eval(parse(text = set))
  mtry_fixed <- sqrt(ncol(data_set))
  for(i in 1:10){
    rf <- randomForest(Diagnosis ~ ., data = data_set, ntree=ntree)
    mod_res <- data.frame(set=set,nr=i,ntree=ntree,oob=rf$err.rate[nrow(rf$err.rate),1])
    tuning_results <- rbind(tuning_results,mod_res)
  }
  
  print(ntree)
}
print(set)
write.csv(tuning_results,paste(set,"_results_rf2.csv",sep=""),row.names = F)


set <- "df4_d"
ntree_seq <- c(seq(280, 299, by = 1)[-1])
tuning_results <- data.frame(set="",nr=0, ntree=0, oob=0)
for (ntree in ntree_seq) {
  data_set <- eval(parse(text = set))
  mtry_fixed <- sqrt(ncol(data_set))
  for(i in 1:10){
    rf <- randomForest(Diagnosis ~ ., data = data_set, ntree=ntree)
    mod_res <- data.frame(set=set,nr=i,ntree=ntree,oob=rf$err.rate[nrow(rf$err.rate),1])
    tuning_results <- rbind(tuning_results,mod_res)
  }
  
  print(ntree)
}
print(set)
write.csv(tuning_results,paste(set,"_results_rf3.csv",sep=""),row.names = F)


set <- "df4_d"
ntree_seq <- c(seq(310, 329, by = 1)[-1])
tuning_results <- data.frame(set="",nr=0, ntree=0, oob=0)
for (ntree in ntree_seq) {
  data_set <- eval(parse(text = set))
  mtry_fixed <- sqrt(ncol(data_set))
  for(i in 1:10){
    rf <- randomForest(Diagnosis ~ ., data = data_set, ntree=ntree)
    mod_res <- data.frame(set=set,nr=i,ntree=ntree,oob=rf$err.rate[nrow(rf$err.rate),1])
    tuning_results <- rbind(tuning_results,mod_res)
  }
  
  print(ntree)
}
print(set)
write.csv(tuning_results,paste(set,"_results_rf4.csv",sep=""),row.names = F)

set <- "df4_d"
data_set <- eval(parse(text = set))
mtry_fixed <- ncol(data_set)

mtry_seq <- c(seq(0, mtry_fixed , by = 10)[-1])
tuning_results <- data.frame(set="",nr=0, mtry=0, oob=0)
for (mtry in mtry_seq) {
  for(i in 1:10){
    rf <- randomForest(Diagnosis ~ ., data = data_set, ntree=680, mtry=mtry)
    mod_res <- data.frame(set=set,nr=i,mtry=mtry,oob=rf$err.rate[nrow(rf$err.rate),1])
    tuning_results <- rbind(tuning_results,mod_res)
  }
  
  print(mtry)
}
print(set)
write.csv(tuning_results,paste(set,"_results_rf_mtry.csv",sep=""),row.names = F)



#There are two ways to find the optimal mtry :
#Apply a similar procedure such that random forest is run 10 times. 
#The optimal number of predictors selected for split is selected for which out of bag error rate stabilizes and reach minimum.
#Experiment with including the (square root of total number of all predictors), 
#(half of this square root value), and (twice of the square root value). 
#And check which mtry returns maximum Area under curve. 
#Thus, for 1000 predictors the number of predictors to select for each node would be 16, 32, and 64 predictors.

###############################
#NTREE
res <- read.csv("df4_a_results_rf2.csv",header = T)
res <- res[-1,]
names <- unique(res$ntree)
# for randomForest
res <- matrix(res$oob,nrow=10)
# for cforest
#res <- matrix(res$kappa,nrow=10)

# Two or more results together
#res2 <- read.csv("df4_a_results_rf2.csv",header = T)
#res2 <- res2[-1,]
#res<- rbind(res,res2)
#res<- res[order(res$ntree,res$nr),]
#res_f <- aggregate(res$oob,by=list(ntree=res$ntree,nr=res$nr),mean)
#res_f<- res_f[order(res_f$ntree,res_f$nr),]
#names <- unique(res_f$ntree)
#res <- res_f
#colnames(res)<-c("ntree","nr","oob")
#res <- matrix(res$oob,nrow=10)

# SD and mean values for iterations
sd_res <- apply(res,2,sd)
mean_res <- apply(res,2,mean)
which(sd_res==min(sd_res))
which(mean_res==min(mean_res))
x_min <-names[which(mean_res==min(mean_res))]

# SD and mean values into dataframe
dat <- data.frame(ntree=names,OOB = mean_res, sd=sd_res)

#Plot version 1
plot(names,mean_res,pch=19,xlab="ntree",ylab="OOB error rate",main="Set A, classifier IV",
     ylim=c(min(mean_res-sd_res),max((mean_res+sd_res))))
lines(rbind(names,names,NA),rbind(mean_res-sd_res,mean_res+sd_res,NA))

#Plot version 2
library(ggplot2)
dat$highlight <- ifelse(dat$ntree %in% c(x_min), "red","black")
ticks <- data.frame (t = c(250, 500, 760, 1000))
png("ntree_tuning_4A.png")
ggplot(dat, aes(x = ntree, y = OOB)) +
  geom_line() +
  geom_ribbon(aes(ymin = OOB - sd,
                  ymax = OOB + sd), alpha = 0.2) +
  labs(y = "OOB error rate")+
  geom_point(aes(x=x_min,y=min(mean_res)), color="red") + 
  scale_x_continuous(breaks=c(ticks$t))
dev.off()

#MTRY
res <- read.csv("df4_a_results_rf_mtry.csv",header = T)
res <- res[-1,]
names <- unique(res$mtry)
res <- matrix(res$oob,nrow=10)
# SD and mean values for iterations
sd_res <- apply(res,2,sd)
mean_res <- apply(res,2,mean)
which(sd_res==min(sd_res))
x_min <-names[which(mean_res==min(mean_res))]
which(mean_res==min(mean_res))
# SD and mean values into dataframe
dat <- data.frame(mtry=names,OOB = mean_res, sd=sd_res)

#Plot version 1
plot(names,mean_res,pch=19,xlab="ntree",ylab="OOB error rate",main="Set A, classifier IV",
     ylim=c(min(mean_res-sd_res),max((mean_res+sd_res))))
lines(rbind(names,names,NA),rbind(mean_res-sd_res,mean_res+sd_res,NA))

#Plot version 2
library(ggplot2)
dat$highlight <- ifelse(dat$ntree %in% c(x_min), "red","black")
ticks <- data.frame (t = c(0, 70, 500,1000,1500))
png("mtry_tuning_4A.png")
ggplot(dat, aes(x = mtry, y = OOB)) +
  geom_line() +
  geom_ribbon(aes(ymin = OOB - sd,
                  ymax = OOB + sd), alpha = 0.2) +
  labs(y = "OOB error rate")+
  geom_point(aes(x=x_min,y=min(mean_res)), color="red") +
  scale_x_continuous(breaks=c(ticks$t))
#scale_x_continuous(breaks=seq(10, 1540, 100)) 
dev.off()