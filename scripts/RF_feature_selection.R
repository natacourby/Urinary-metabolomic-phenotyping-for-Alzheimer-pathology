##################################################################################################
# Feature selection methods - Random Forests
# 1) Choose set and classifier
# Set A - Metabolites only, set B - Metabolites and SNPs, set C - Metabolites, SNPs and covariates, set D - Metabolites and covariates
# Classifier I - AD/CTL/cMCI/sMCI, Classifier II - AD+cMCI/CTL+sMCI, Classifier III - AD/CTL
# 2) Tune RF parameters for chosen dataset (set and classifier)
# 3) Run final RF model with importance parameter to rank metabolites
##################################################################################################

library(randomForest)
library(dplyr)
library(caret)
library(Boruta)
library(ggplot2)
source("./scripts/feature_selection.R")

nr <- 500

feature_selection_RF_model <- function(set, classifier){
  df_selected <- feature_selection_data_preparation(set,classifier)
  err_rates <- c(0)
  for(i in 1:10){
    rf=randomForest(Diagnosis ~ ., data = df_selected)
    err_rates <- rbind(err_rates,rf$err.rate[,1])
  }
  err_rates <- err_rates[-1,]
  result <- apply(err_rates,2,mean)
  return (result)
}

# 1) Choose set and classifier 
# Set A: Metabolites only 
mean_res1_a <- feature_selection_RF_model("A","I")
mean_res2_a <- feature_selection_RF_model("A","II")
mean_res4_a <- feature_selection_RF_model("A","III")
# Set B: Metabolites and SNPs
mean_res1_b <- feature_selection_RF_model("B","I")
mean_res2_b <- feature_selection_RF_model("B","II")
mean_res4_b <- feature_selection_RF_model("B","III")
# Set C: Metabolites, SNPs and covariates
mean_res1_c <- feature_selection_RF_model("C","I")
mean_res2_c <- feature_selection_RF_model("C","II")
mean_res4_c <- feature_selection_RF_model("C","III")
# Set D: Metabolites and covariates
mean_res1_d <- feature_selection_RF_model("D","I")
mean_res2_d <- feature_selection_RF_model("D","II")
mean_res4_d <- feature_selection_RF_model("D","III")

err_rates <- cbind(mean_res1_a,mean_res2_a)
err_rates <- cbind(err_rates,mean_res4_a)
err_rates <- cbind(err_rates,mean_res1_b)
err_rates <- cbind(err_rates,mean_res2_b)
err_rates <- cbind(err_rates,mean_res4_b)
err_rates <- cbind(err_rates,mean_res1_c)
err_rates <- cbind(err_rates,mean_res2_c)
err_rates <- cbind(err_rates,mean_res4_c)
err_rates <- cbind(err_rates,mean_res1_d)
err_rates <- cbind(err_rates,mean_res2_d)
err_rates <- cbind(err_rates,mean_res4_d)

# Colours for sets (A, B, C, D)
cols <- c("#f58231",
          "#f58231",
          "#f58231",
          "#000080",
          "#000080",
          "#000080",
          "#808000",
          "#808000",
          "#808000",
          "#e6194b",
          "#e6194b",
          "#e6194b")

# Line types for classifiers (I, II, III)
ltys <- c(2,3,1,2,3,1,2,3,1,2,3,1)

png("./images/RF_all_sets.png")
matplot(1:nr , err_rates, col=cols,type="l",lty=ltys, ylab="OOB error",xlab="trees")
legend("topright",
       legend=c("Set A, Classifier I","Set A, Classifier II","Set A, Classifier III",
                "Set B, Classifier I","Set B, Classifier II","Set B, Classifier III",
                "Set C, Classifier I","Set C, Classifier II","Set C, Classifier III",
                "Set D, Classifier I","Set D, Classifier II","Set D, Classifier III")
       ,lty=ltys, 
       col=cols,cex=0.8)
dev.off()

min(err_rates)
min(mean_res4_a)
min(mean_res4_d)
# Set D, Classifier III - smallest OOB = 0.02847458
which(mean_res4_d==min(mean_res4_d))
#477

# 2) Tune RF parameters
data_set <- feature_selection_data_preparation("D","III")

# NTREE
tuning_results <- data.frame(set="",nr=0, ntree=0, oob=0)
control <- trainControl(method="oob")
ntree_seq <- c(seq(0, 1000, by = 10)[-1])

tuning_results_ntree <- data.frame(set="",nr=0, ntree=0, oob=0)
for (ntree in ntree_seq) {
  mtry_fixed <- sqrt(ncol(data_set))
  for(i in 1:10){
    rf <- randomForest(Diagnosis ~ ., data = data_set, ntree=ntree)
    mod_res <- data.frame(set=set,nr=i,ntree=ntree,oob=rf$err.rate[nrow(rf$err.rate),1])
    tuning_results_ntree <- rbind(tuning_results_ntree,mod_res)
  }
  
  print(ntree)
}
write.csv(tuning_results_ntree,"./results/RF_SetD_ClassifierIII_results_ntree.csv",row.names = F)

res <- matrix(tuning_results_ntree$oob,nrow=10)
mean_res <- apply(res,2,mean)
i <- which(mean_res==min(mean_res))
ntree_selected <-i*10

# MTRY 
mtry_fixed <- ncol(data_set)
mtry_seq <- c(seq(0, mtry_fixed , by = 10)[-1])
tuning_results_mtry <- data.frame(set="",nr=0, mtry=0, oob=0)
for (mtry in mtry_seq) {
  for(i in 1:10){
    rf <- randomForest(Diagnosis ~ ., data = data_set, ntree=ntree_selected, mtry=mtry)
    mod_res <- data.frame(set=set,nr=i,mtry=mtry,oob=rf$err.rate[nrow(rf$err.rate),1])
    tuning_results_mtry <- rbind(tuning_results_mtry,mod_res)
  }
  print(mtry)
}
write.csv(tuning_results_mtry,"./results/RF_SetD_ClassifierIII_results_mtry.csv",row.names = F)

res <- matrix(tuning_results_mtry$oob,nrow=10)
mean_res <- apply(res,2,mean)
i <- which(mean_res==min(mean_res))
mtry_selected <-i*10

#NTREE
res <- tuning_results_ntree
res <- res[-1,]
names <- unique(res$ntree)
# for randomForest
res <- matrix(res$oob,nrow=10)

# SD and mean values for iterations
sd_res <- apply(res,2,sd)
mean_res <- apply(res,2,mean)
which(sd_res==min(sd_res))
which(mean_res==min(mean_res))
x_min <-names[which(mean_res==min(mean_res))]

# SD and mean values into dataframe
dat <- data.frame(ntree=names,OOB = mean_res, sd=sd_res)

#Plot version 1
plot(names,mean_res,pch=19,xlab="ntree",ylab="OOB error rate",main="Set D, ClassifierIII",
     ylim=c(min(mean_res-sd_res),max((mean_res+sd_res))))
lines(rbind(names,names,NA),rbind(mean_res-sd_res,mean_res+sd_res,NA))

#Plot version 2
dat$highlight <- ifelse(dat$ntree %in% c(x_min), "red","black")
ticks <- data.frame (t = c(250, 500, x_min, 1000))
png("./images/RF_ntree_tuning.png")
ggplot(dat, aes(x = ntree, y = OOB)) +
  geom_line() +
  geom_ribbon(aes(ymin = OOB - sd,
                  ymax = OOB + sd), alpha = 0.2) +
  labs(y = "OOB error")+
  geom_point(aes(x=x_min,y=min(mean_res)), color="red") + 
  scale_x_continuous(breaks=c(ticks$t)) + 
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=16))
dev.off()

#MTRY
res <- tuning_results_mtry
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
plot(names,mean_res,pch=19,xlab="ntree",ylab="OOB error",main="Set D, ClassifierIII",
     ylim=c(min(mean_res-sd_res),max((mean_res+sd_res))))
lines(rbind(names,names,NA),rbind(mean_res-sd_res,mean_res+sd_res,NA))

#Plot version 2
dat$highlight <- ifelse(dat$mtry %in% c(x_min), "red","black")
ticks <- data.frame (t = c(0, x_min, 500,1000,1500))
png("./images/RF_mtree_tuning.png")
ggplot(dat, aes(x = mtry, y = OOB)) +
  geom_line() +
  geom_ribbon(aes(ymin = OOB - sd,
                  ymax = OOB + sd), alpha = 0.2) +
  labs(y = "OOB error")+
  geom_point(aes(x=x_min,y=min(mean_res)), color="red") +
  scale_x_continuous(breaks=c(ticks$t)) + 
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=16))
dev.off()

# 3) Run final RF model with importance parameter to rank metabolites
model <- randomForest(Diagnosis ~ ., data = data_set, ntree=ntree_selected, mtry=mtry_selected, importance=TRUE)
importance    <- randomForest::importance(model)
varImportance <- data.frame(Variables = row.names(importance), Importance = round(importance[ ,'MeanDecreaseGini'],2))

#Create a rank variable based on importance
rankImportance <- varImportance %>% mutate(Rank = paste0('#',dense_rank(desc(Importance))))

# Boruta algorithm
boruta.train <- Boruta(Diagnosis~., data = data_set, doTrace = 2)

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

write.csv(getSelectedAttributes(boruta.train, withTentative = T),"./results/RF_Boruta_features_ranked.txt", row.names = FALSE)
write.csv(varImportance,"./results/RF_features_ranked.txt", row.names = FALSE)

