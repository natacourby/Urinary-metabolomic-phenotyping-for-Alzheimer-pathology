################################################################################
# RF tuning parameters ntree and mtry
################################################################################
set <- "df4_d"
library(randomForest)
source("/hps/nobackup/ma/natalja/data/mQTL_data/scripts/load_data.R")

# NTREE
tuning_results <- data.frame(set="",nr=0, ntree=0, oob=0)
control <- trainControl(method="oob")
ntree_seq <- c(seq(0, 1000, by = 10)[-1])

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
write.csv(tuning_results,paste(set,"_results_rf.csv",sep=""),row.names = F)

res <- matrix(tuning_results$oob,nrow=10)
mean_res <- apply(res,2,mean)
ntree_selected <-names[which(mean_res==min(mean_res))]

# MTRY 
data_set <- eval(parse(text = set))
mtry_fixed <- ncol(data_set)
mtry_seq <- c(seq(0, mtry_fixed , by = 10)[-1])
tuning_results <- data.frame(set="",nr=0, mtry=0, oob=0)
for (mtry in mtry_seq) {
  for(i in 1:10){
    rf <- randomForest(Diagnosis ~ ., data = data_set, ntree=ntree_selected, mtry=mtry)
    mod_res <- data.frame(set=set,nr=i,mtry=mtry,oob=rf$err.rate[nrow(rf$err.rate),1])
    tuning_results <- rbind(tuning_results,mod_res)
  }
  print(mtry)
}
write.csv(tuning_results,paste(set,"_results_rf_mtry.csv",sep=""),row.names = F)

