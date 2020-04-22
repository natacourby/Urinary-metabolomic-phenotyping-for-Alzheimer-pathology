# Script plots RF parameters ntree and mtry
set <- "4_d"
set_name <- "Set D, classifier IV"
#NTREE
res <- read.csv(paste("df",set,"_results_rf.csv",sep=""),header = T)
res <- res[-1,]
names <- unique(res$ntree)
# for randomForest
res <- matrix(res$oob,nrow=10)
# for cforest
#res <- matrix(res$kappa,nrow=10)


# SD and mean values for iterations
sd_res <- apply(res,2,sd)
mean_res <- apply(res,2,mean)
which(sd_res==min(sd_res))
which(mean_res==min(mean_res))
x_min <-names[which(mean_res==min(mean_res))]

# SD and mean values into dataframe
dat <- data.frame(ntree=names,OOB = mean_res, sd=sd_res)

#Plot version 1
plot(names,mean_res,pch=19,xlab="ntree",ylab="OOB error rate",main=set_name,
     ylim=c(min(mean_res-sd_res),max((mean_res+sd_res))))
lines(rbind(names,names,NA),rbind(mean_res-sd_res,mean_res+sd_res,NA))

#Plot version 2
library(ggplot2)
dat$highlight <- ifelse(dat$ntree %in% c(x_min), "red","black")
ticks <- data.frame (t = c(250, 500, x_min, 1000))
png(paste("ntree_tuning_",set,".png",sep=""))
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
res <- read.csv(paste("df",set,"_results_rf_mtry.csv",sep=""),header = T)
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
plot(names,mean_res,pch=19,xlab="ntree",ylab="OOB error",main=set_name,
     ylim=c(min(mean_res-sd_res),max((mean_res+sd_res))))
lines(rbind(names,names,NA),rbind(mean_res-sd_res,mean_res+sd_res,NA))

#Plot version 2
library(ggplot2)
dat$highlight <- ifelse(dat$ntree %in% c(x_min), "red","black")
ticks <- data.frame (t = c(0, x_min, 500,1000,1500))
png(paste("mtry_tuning_",set,".png",sep=""))
ggplot(dat, aes(x = mtry, y = OOB)) +
  geom_line() +
  geom_ribbon(aes(ymin = OOB - sd,
                  ymax = OOB + sd), alpha = 0.2) +
  labs(y = "OOB error")+
  geom_point(aes(x=x_min,y=min(mean_res)), color="red") +
  scale_x_continuous(breaks=c(ticks$t)) + 
  theme(axis.text=element_text(size=14),
         axis.title=element_text(size=16))
  #scale_x_continuous(breaks=seq(10, 1540, 100)) 
dev.off()