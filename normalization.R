library(sva)
library(ggfortify)
source("/home/natalja/Documents/R_projects/EigenMS.R")
QN_normalize <- function(file_name, output_file_name){
  met <-read.table(file_name,header=T,sep="\t")
  for (i in 1:nrow(met)){
    mat <- met[i,]
    mat = t(apply(mat, 1, rank, ties.method = "average"));
    mat = qnorm(mat / (ncol(met)+1))
    met[i,] <- mat 
  }
  write.table(met,output_file_name,sep="\t")
}
################################################################################
# LC-MS SERUM
# 1) Remove biases (cohort and centre) by ComBat
# 2) EigenMS
# 3) QN
################################################################################
bases <- c("SHPOS","SLPOS","SLNEG")
for(base in bases){
  main_dir <- "/home/natalja/Documents/EMIF/KCL/NPC/NPC/original_data/"

  file_name <- paste(main_dir,base,"_original.csv",sep="")
  output_file_name <- paste("/home/natalja/Documents/EMIF/KCL/NPC/NPC/normalized_data/EigenMS_removed_cohort_centre/",base,"_EigenMS.tsv",sep="")
  qn_output_file_name <- paste("/home/natalja/Documents/EMIF/KCL/NPC/NPC/normalized_data/EigenMS_removed_cohort_centre/",base,"_EigenMS_QN.tsv",sep="")
  
  met <-read.table(file_name,header=T,sep=",")
  rownames(met) <- met$name
  met <- met[,-(1:3)]
  met <- met[complete.cases(met), ]
  if (base =="SHPOS"){
    qn_met <- met[,substr(colnames(met),1,2)=="QC"]
    met <- met[,substr(colnames(met),1,2)!="QC"]
  } else {
    qn_met <- met[,substr(colnames(met),1,2)=="SR"]
    met <- met[,substr(colnames(met),1,2)!="SR"]
  }
  
  cov_file_name <- "/home/natalja/Documents/EMIF/KCL/Covariates_all_vertical_text2.txt"
 # cov_file_name <- "/home/natalja/Documents/EMIF/KCL/Covariates_all_categories.txt"
  
  cov_file_name_gm <- "/home/natalja/Documents/EMIF/KCL/mQTL/mQTL_data/Covariates/Covariates_all_gm_vertical.txt"
  
  cov <- read.table(cov_file_name,header = T,sep="\t", stringsAsFactors = FALSE)
  cov_gm <- read.table(cov_file_name_gm,header = T,sep="\t", stringsAsFactors = FALSE)

  rownames(cov_gm) <- cov_gm$Id
  rownames(cov) <- cov$Sample
  
  #cov<-cov[rownames(cov) %in% rownames(cov_gm),]
  
  
  rownames(cov) <- cov$Sample
  cov<-cov[colnames(met),]
  
  e <- met[,names(met) %in% rownames(cov)]
  cov <- cov[rownames(cov) %in% colnames(e),]
  cov <- cov[colnames(e),]
  
  
  
  
  
  

  ddata <- as.matrix(e)
###################################  
# Original PCA
###################################    
  scaleShift=abs(min(t(ddata), na.rm = TRUE))+1 
  df = log(t(ddata)+scaleShift)
  cdata <- as.data.frame(df)
  cdata$Diagnosis <- cov$Diagnosis
  cdata$Cohort <- cov$Cohort
  cdata$Centre <- cov$Centre
  cdata$Gender <- cov$Gender
  
  setwd("/home/natalja/Documents/EMIF/images_mQTL/PCA/")
  png(paste("PCA_",base,"_Diagnosis_original.png",sep=""), width=950, height=400)
  autoplot(prcomp(df), data = cdata, colour = 'Diagnosis', main=base)
  dev.off()
  png(paste("PCA_",base,"_Cohort_original.png",sep=""), width=950, height=400)
  autoplot(prcomp(df), data = cdata, colour = 'Cohort', main=base)
  dev.off()
  png(paste("PCA_",base,"_Centre_original.png",sep=""), width=950, height=400)
  autoplot(prcomp(df), data = cdata, colour = 'Centre', main=base)
  dev.off()
  
###################################  
# Remove biases with ComBat
###################################    
  
  # Remove cohort bias
  batch = as.factor(as.character(cov$BIAS_PK))
  modcombat = model.matrix(~1, data=cov)
  combat_data = ComBat(dat=ddata, batch=batch, mod=modcombat, par.prior=TRUE, prior.plots=FALSE)
  ddata <- combat_data
  # Remove centre bias
  batch = as.factor(as.character(cov$Centre))
  modcombat = model.matrix(~1, data=cov)
  ddata = ComBat(dat=ddata, batch=batch, mod=modcombat, par.prior=TRUE, prior.plots=FALSE)

  X <- model.matrix(~1 + Centre,data = cov)
  output <- lm(t(original_ddata) ~ X)
  lr_ddata <- t(output$residuals)
  
  
###################################  
# Removed biases PCA
###################################    
  scaleShift=abs(min(t(lr_ddata), na.rm = TRUE))+1
  df = log(t(lr_ddata)+scaleShift)
  cdata <- as.data.frame(df)
  cdata$Diagnosis <- cov$Diagnosis
  cdata$Cohort <- cov$Cohort
  cdata$Centre <- cov$Centre
  cdata$Gender <- cov$Gender
  png(paste("PCA_",base,"_Diagnosis_removed_cohort_centre.png",sep=""), width=950, height=400)
  autoplot(prcomp(df), data = cdata, colour = 'Diagnosis', main=base)
  dev.off()
  #png(paste("PCA_",base,"_Cohort_removed_cohort_centre.png",sep=""), width=950, height=400)
  #autoplot(prcomp(df), data = cdata, colour = 'Cohort', main=base)
  #dev.off()
  png(paste("PCA_",base,"_Centre_removed_cohort_centre.png",sep=""), width=950, height=400)
  autoplot(prcomp(df), data = cdata, colour = 'Centre', main=base)
  dev.off()

###################################  
# EigenMS
###################################    
  
  #grps <- data.frame(Diagnosis=as.factor(cov$EigenMS_factor_to.preserve),Centre=as.numeric(as.factor(cov$Centre)))
  grps <- data.frame(diagnosis = as.factor((as.numeric(cov$EigenMS_factor_to.preserve))),centre = as.factor((as.numeric(as.factor(cov$Centre)))) )
  rownames(grps) <- cov$Sample
  
  grps <- as.factor(cov$EigenMS_factor_to.preserve)
  names(grps) <- cov$Sample

  scaleShift=abs(min(original_ddata, na.rm = TRUE))+1  
  m_logInts = log(original_ddata+scaleShift)

  m_nummiss = sum(is.na(m_logInts)) 

  m_numtot = dim(m_logInts)[1] * dim(m_logInts)[2]
  m_percmiss = m_nummiss/m_numtot  
  m_prot.info = cbind(rownames(ddata),rownames(ddata)) 
  
  m_ints_eig1 = eig_norm1(m=m_logInts,treatment=grps,prot.info=m_prot.info,write_to_file = output_file_name)
  #m_ints_eig1$h.c = 4
  m_ints_norm1 = eig_norm2(rv=m_ints_eig1) 
  
  eigenMS2_ddata <- m_ints_norm1$norm_m  
  write.table(ddata,output_file_name,sep="\t")
  
###################################  
# EigenMS PCA
################################### 
  #ddata <- ddata[,colnames(ddata) %in% rownames(cov)]

  scaleShift=abs(min(t(eigenMS1_ddata), na.rm = TRUE))+1 
  df = log(t(eigenMS1_ddata)+scaleShift)
  cdata <- as.data.frame(df)
  cdata$Diagnosis <- cov$Diagnosis
  cdata$Cohort <- cov$Cohort
  cdata$Centre <- cov$Centre
  cdata$Gender <- cov$Gender
  png(paste("PCA_",base,"_Diagnosis_EigenMS.png",sep=""), width=950, height=400)
  autoplot(prcomp(df), data = cdata, colour = 'Diagnosis', main=base)
  dev.off()
  #png(paste("PCA_",base,"_Cohort_EigenMS.png",sep=""), width=950, height=400)
  #autoplot(prcomp(df), data = cdata, colour = 'Cohort', main=base)
  #dev.off()
  png(paste("PCA_",base,"_Centre_EigenMS.png",sep=""), width=950, height=400)
  autoplot(prcomp(df), data = cdata, colour = 'Centre', main=base)
  dev.off()
  
###################################  
# QN
###################################   

  # Run and save QN results
  QN_normalize(output_file_name,qn_output_file_name)
  
###################################  
# QN PCA
################################### 

  met <-read.table(qn_output_file_name,header=T,sep="\t")
  cov_file_name <- "/home/natalja/Documents/EMIF/KCL/Covariates_all_vertical_text2.txt"
  cov <- read.table(cov_file_name,header = T,sep="\t", stringsAsFactors = FALSE)
  rownames(cov) <- cov$Sample
  cov<-cov[colnames(met),]
  
  e <- met[,names(met) %in% rownames(cov)]
  cov <- cov[rownames(cov) %in% colnames(e),]
  cov <- cov[colnames(e),]
  
  ddata <- as.matrix(e)

  scaleShift=abs(min(t(ddata), na.rm = TRUE))+1 
  df = log(t(ddata)+scaleShift)
  cdata <- as.data.frame(df)
  cdata$Diagnosis <- cov$Diagnosis
  cdata$Cohort <- cov$Cohort
  cdata$Centre <- cov$Centre
  cdata$Gender <- cov$Gender
  png(paste("PCA_",base,"_Diagnosis_QN.png",sep=""), width=950, height=400)
  autoplot(prcomp(df), data = cdata, colour = 'Diagnosis', main=base)
  dev.off()
  #png(paste("PCA_",base,"_Cohort_QN.png",sep=""), width=950, height=400)
  #autoplot(prcomp(df), data = cdata, colour = 'Cohort', main=base)
  #dev.off()
  png(paste("PCA_",base,"_Centre_QN.png",sep=""), width=950, height=400)
  autoplot(prcomp(df), data = cdata, colour = 'Centre', main=base)
  dev.off()
  
}  

################################################################################
# LC-MS URINE
# 1) EigenMS 
# 2) QN 
################################################################################
bases <- c("UHPOS","URPOS","URNEG")
for(base in bases){
  main_dir <- "/home/natalja/Documents/EMIF/KCL/NPC/NPC/original_data/"
  
  file_name <- paste(main_dir,base,"_original.csv",sep="")
  output_file_name <- paste("/home/natalja/Documents/EMIF/KCL/NPC/NPC/normalized_data/EigenMS_removed_cohort_centre/",base,"_EigenMS.tsv",sep="")
  qn_output_file_name <- paste("/home/natalja/Documents/EMIF/KCL/NPC/NPC/normalized_data/EigenMS_removed_cohort_centre/",base,"_EigenMS_QN.tsv",sep="")
  
  met <-read.table(file_name,header=T,sep=",")
  rownames(met) <- met$name
  met <- met[,-(1:3)]
  met <- met[complete.cases(met), ]
  qn_met <- met[,substr(colnames(met),1,2)=="SR"]
  met <- met[,substr(colnames(met),1,2)!="SR"]

  cov_file_name <- "/home/natalja/Documents/EMIF/KCL/Covariates_all_vertical_text.txt"
  cov <- read.table(cov_file_name,header = T,sep="\t", stringsAsFactors = FALSE)
  rownames(cov) <- cov$Sample
  cov<-cov[colnames(met),]
  
  e <- met[,names(met) %in% rownames(cov)]
  cov <- cov[rownames(cov) %in% colnames(e),]
  cov <- cov[colnames(e),]
  
  ddata <- as.matrix(e)
  
  ###################################  
  # Original PCA
  ###################################    
  scaleShift=abs(min(t(ddata), na.rm = TRUE))+1 
  df = log(t(ddata)+scaleShift)
  cdata <- as.data.frame(df)
  cdata$Diagnosis <- cov$Diagnosis
  cdata$Cohort <- cov$Cohort
  cdata$Centre <- cov$Centre
  cdata$Gender <- cov$Gender
  
  setwd("/home/natalja/Documents/EMIF/images_mQTL/")
  png(paste("PCA_",base,"_Diagnosis_original.png",sep=""), width=950, height=400)
  autoplot(prcomp(df), data = cdata, colour = 'Diagnosis', main=base)
  dev.off()
  png(paste("PCA_",base,"_Cohort_original.png",sep=""), width=950, height=400)
  autoplot(prcomp(df), data = cdata, colour = 'Cohort', main=base)
  dev.off()
  png(paste("PCA_",base,"_Centre_original.png",sep=""), width=950, height=400)
  autoplot(prcomp(df), data = cdata, colour = 'Centre', main=base)
  dev.off()
  
  ###################################  
  # EigenMS
  ###################################    
  
  grps <- as.factor(cov$EigenMS_factor)
  
  scaleShift=abs(min(ddata, na.rm = TRUE))+1  
  m_logInts = log(ddata+scaleShift)
  
  m_nummiss = sum(is.na(m_logInts)) 
  
  m_numtot = dim(m_logInts)[1] * dim(m_logInts)[2]
  m_percmiss = m_nummiss/m_numtot  
  m_prot.info = cbind(rownames(ddata),rownames(ddata)) 
  
  m_ints_eig1 = eig_norm1(m=m_logInts,treatment=grps,prot.info=m_prot.info,write_to_file = output_file_name)
  m_ints_norm1 = eig_norm2(rv=m_ints_eig1) 
  
  ddata <- m_ints_norm1$norm_m  
  write.table(ddata,output_file_name,sep="\t")
  
  ###################################  
  # EigenMS PCA
  ################################### 
  ddata <- ddata[,colnames(ddata) %in% rownames(cov)]
  
  scaleShift=abs(min(t(ddata), na.rm = TRUE))+1 
  df = log(t(ddata)+scaleShift)
  cdata <- as.data.frame(df)
  cdata$Diagnosis <- cov$Diagnosis
  cdata$Cohort <- cov$Cohort
  cdata$Centre <- cov$Centre
  cdata$Gender <- cov$Gender
  png(paste("PCA_",base,"_Diagnosis_EigenMS.png",sep=""), width=950, height=400)
  autoplot(prcomp(df), data = cdata, colour = 'Diagnosis', main=base)
  dev.off()
  png(paste("PCA_",base,"_Cohort_EigenMS.png",sep=""), width=950, height=400)
  autoplot(prcomp(df), data = cdata, colour = 'Cohort', main=base)
  dev.off()
  png(paste("PCA_",base,"_Centre_EigenMS.png",sep=""), width=950, height=400)
  autoplot(prcomp(df), data = cdata, colour = 'Centre', main=base)
  dev.off()
  
  ###################################  
  # QN
  ###################################   
  
  # Run and save QN results
  QN_normalize(output_file_name,qn_output_file_name)
  
  ###################################  
  # QN PCA
  ################################### 
  
  met <-read.table(qn_output_file_name,header=T,sep="\t")
  cov_file_name <- "/home/natalja/Documents/EMIF/KCL/Covariates_all_vertical_text.txt"
  cov <- read.table(cov_file_name,header = T,sep="\t", stringsAsFactors = FALSE)
  rownames(cov) <- cov$Sample
  cov<-cov[colnames(met),]
  
  e <- met[,names(met) %in% rownames(cov)]
  cov <- cov[rownames(cov) %in% colnames(e),]
  cov <- cov[colnames(e),]
  
  ddata <- as.matrix(e)
  
  scaleShift=abs(min(t(ddata), na.rm = TRUE))+1 
  df = log(t(ddata)+scaleShift)
  cdata <- as.data.frame(df)
  cdata$Diagnosis <- cov$Diagnosis
  cdata$Cohort <- cov$Cohort
  cdata$Centre <- cov$Centre
  cdata$Gender <- cov$Gender
  png(paste("PCA_",base,"_Diagnosis_QN.png",sep=""), width=950, height=400)
  autoplot(prcomp(df), data = cdata, colour = 'Diagnosis', main=base)
  dev.off()
  png(paste("PCA_",base,"_Cohort_QN.png",sep=""), width=950, height=400)
  autoplot(prcomp(df), data = cdata, colour = 'Cohort', main=base)
  dev.off()
  png(paste("PCA_",base,"_Centre_QN.png",sep=""), width=950, height=400)
  autoplot(prcomp(df), data = cdata, colour = 'Centre', main=base)
  dev.off()
  
} 
################################################################################
# NMR SERUM
# 1) Remove biases (cohort and centre) by ComBat
# 2) QN
################################################################################
bases <- c("NMR_CPMG_Serum","NMR_Nosey_Serum","NMR_BI_LISA") 
for(base in bases){
  main_dir <- "/home/natalja/Documents/EMIF/KCL/NPC/NPC/NMR/original_data/"
  
  file_name <- paste(main_dir,base,"_original.txt",sep="")
  output_file_name <- paste("/home/natalja/Documents/EMIF/KCL/NPC/NPC/NMR/normalized_data/EigenMS_removed_cohort_centre/",base,"_rb.tsv",sep="")
  qn_output_file_name <- paste("/home/natalja/Documents/EMIF/KCL/NPC/NPC/NMR/normalized_data/EigenMS_removed_cohort_centre/",base,"_rb_QN.tsv",sep="")
  
  met <-read.table(file_name,header=T,sep="\t",stringsAsFactors = F)
  if (base =="NMR_BI_LISA"){
    dup <- met[met$Sampling.ID %in% met[duplicated(met$Sampling.ID),]$Sampling.ID,]
    met <- met[!met$Sampling.ID %in% met[duplicated(met$Sampling.ID),]$Sampling.ID,] 

    library(plyr)
    met <- rbind(met,ddply(dup,"Sampling.ID",numcolwise(mean)))
    
    rownames(met) <- met$Sampling.ID
    
    met <- met[,-1]
    met <- t(met)
  }
     

  cov_file_name <- "/home/natalja/Documents/EMIF/KCL/Covariates_all_vertical_text2.txt"
  cov <- read.table(cov_file_name,header = T,sep="\t", stringsAsFactors = FALSE)
  rownames(cov) <- cov$Sample
  cov<-cov[colnames(met),]
  
  e <- met[,colnames(met) %in% rownames(cov)]
  cov <- cov[rownames(cov) %in% colnames(e),]
  cov <- cov[colnames(e),]
  
  ddata <- as.matrix(e)
  
  ###################################  
  # Original PCA
  ###################################    
  scaleShift=abs(min(t(ddata), na.rm = TRUE))+1 
  df = log(t(ddata)+scaleShift)
  cdata <- as.data.frame(df)
  cdata$Diagnosis <- cov$Diagnosis
  cdata$Cohort <- cov$Cohort
  cdata$Centre <- cov$Centre
  cdata$Gender <- cov$Gender
  
  setwd("/home/natalja/Documents/EMIF/images_mQTL/PCA/")
  png(paste("PCA_",base,"_Diagnosis_original.png",sep=""), width=950, height=400)
  autoplot(prcomp(df), data = cdata, colour = 'Diagnosis', main=base)
  dev.off()
  png(paste("PCA_",base,"_Cohort_original.png",sep=""), width=950, height=400)
  autoplot(prcomp(df), data = cdata, colour = 'Cohort', main=base)
  dev.off()
  png(paste("PCA_",base,"_Centre_original.png",sep=""), width=950, height=400)
  autoplot(prcomp(df), data = cdata, colour = 'Centre', main=base)
  dev.off()
  
  ###################################  
  # Remove biases with ComBat
  ###################################    
  
  # Remove cohort bias
  batch = cov$Cohort
  modcombat = model.matrix(~1, data=cov)
  combat_data = ComBat(dat=ddata, batch=batch, mod=modcombat, par.prior=TRUE, prior.plots=FALSE)
  
  write.table(combat_data,paste("/home/natalja/Documents/EMIF/KCL/NPC/NPC/NMR/original_data/",base,"_cohort_bias_removed.txt",sep=""),sep="\t",row.names = TRUE)
  
  ddata <- combat_data
  
  # Remove centre bias
  batch = cov$Centre
  modcombat = model.matrix(~1, data=cov)
  combat_data2 = ComBat(dat=ddata, batch=batch, mod=modcombat, par.prior=TRUE, prior.plots=FALSE)
  
  write.table(combat_data2,paste("/home/natalja/Documents/EMIF/KCL/NPC/NPC/NMR/original_data/",base,"_cohort_centre_bias_removed.txt",sep=""),sep="\t",row.names = TRUE)
  
  ddata <- combat_data2
  
  ###################################  
  # Removed biases PCA
  ###################################    
  scaleShift=abs(min(t(ddata), na.rm = TRUE))+1 
  df = log(t(ddata)+scaleShift)
  cdata <- as.data.frame(df)
  cdata$Diagnosis <- cov$Diagnosis
  cdata$Cohort <- cov$Cohort
  cdata$Centre <- cov$Centre
  cdata$Gender <- cov$Gender
  
  setwd("/home/natalja/Documents/EMIF/images_mQTL/PCA/")
  png(paste("PCA_",base,"_Diagnosis_removed_cohort_centre.png",sep=""), width=950, height=400)
  autoplot(prcomp(df), data = cdata, colour = 'Diagnosis', main=base)
  dev.off()
  png(paste("PCA_",base,"_Cohort_removed_cohort_centre.png",sep=""), width=950, height=400)
  autoplot(prcomp(df), data = cdata, colour = 'Cohort', main=base)
  dev.off()
  png(paste("PCA_",base,"_Centre_removed_cohort_centre.png",sep=""), width=950, height=400)
  autoplot(prcomp(df), data = cdata, colour = 'Centre', main=base)
  dev.off()
  
  write.table(ddata,output_file_name,sep="\t")
  
  ###################################  
  # QN
  ################################### 
  
  QN_normalize(output_file_name,qn_output_file_name)
  
  ###################################  
  # QN PCA
  ################################### 
  
  met <-read.table(qn_output_file_name,header=T,sep="\t")
  cov_file_name <- "/home/natalja/Documents/EMIF/KCL/Covariates_all_vertical_text2.txt"
  cov <- read.table(cov_file_name,header = T,sep="\t", stringsAsFactors = FALSE)
  rownames(cov) <- cov$Sample
  cov<-cov[colnames(met),]
  
  e <- met[,names(met) %in% rownames(cov)]
  cov <- cov[rownames(cov) %in% colnames(e),]
  cov <- cov[colnames(e),]
  
  ddata <- as.matrix(e)
  
  scaleShift=abs(min(t(ddata), na.rm = TRUE))+1 
  df = log(t(ddata)+scaleShift)
  cdata <- as.data.frame(df)
  cdata$Diagnosis <- cov$Diagnosis
  cdata$Cohort <- cov$Cohort
  cdata$Centre <- cov$Centre
  cdata$Gender <- cov$Gender
  png(paste("PCA_",base,"_Diagnosis_QN.png",sep=""), width=950, height=400)
  autoplot(prcomp(df), data = cdata, colour = 'Diagnosis', main=base)
  dev.off()
  png(paste("PCA_",base,"_Cohort_QN.png",sep=""), width=950, height=400)
  autoplot(prcomp(df), data = cdata, colour = 'Cohort', main=base)
  dev.off()
  png(paste("PCA_",base,"_Centre_QN.png",sep=""), width=950, height=400)
  autoplot(prcomp(df), data = cdata, colour = 'Centre', main=base)
  dev.off()
}  
################################################################################
# NMR URINE
# 1) QN
################################################################################
bases <- c("NMR_Nosey_Urine") 
for(base in bases){
  main_dir <- "/home/natalja/Documents/EMIF/KCL/NPC/NPC/NMR/original_data/"
  
  file_name <- paste(main_dir,base,"_original.txt",sep="")
  output_file_name <- paste("/home/natalja/Documents/EMIF/KCL/NPC/NPC/NMR/normalized_data/EigenMS_removed_cohort_centre/",base,"_EigenMS.tsv",sep="")
  qn_output_file_name <- paste("/home/natalja/Documents/EMIF/KCL/NPC/NPC/NMR/normalized_data/EigenMS_removed_cohort_centre/",base,"_EigenMS_QN.tsv",sep="")
  
  met <-read.table(file_name,header=T,sep="\t")
  
  cov_file_name <- "/home/natalja/Documents/EMIF/KCL/Covariates_all_vertical_text.txt"
  cov <- read.table(cov_file_name,header = T,sep="\t", stringsAsFactors = FALSE)
  rownames(cov) <- cov$Sample
  cov<-cov[colnames(met),]
  
  e <- met[,names(met) %in% rownames(cov)]
  cov <- cov[rownames(cov) %in% colnames(e),]
  cov <- cov[colnames(e),]
  
  ddata <- as.matrix(e)
  
  ###################################  
  # Original PCA
  ###################################    
  scaleShift=abs(min(t(ddata), na.rm = TRUE))+1 
  df = log(t(ddata)+scaleShift)
  cdata <- as.data.frame(df)
  cdata$Diagnosis <- cov$Diagnosis
  cdata$Cohort <- cov$Cohort
  cdata$Centre <- cov$Centre
  cdata$Gender <- cov$Gender
  
  setwd("/home/natalja/Documents/EMIF/images_mQTL/")
  png(paste("PCA_",base,"_Diagnosis_original.png",sep=""), width=950, height=400)
  autoplot(prcomp(df), data = cdata, colour = 'Diagnosis', main=base)
  dev.off()
  png(paste("PCA_",base,"_Cohort_original.png",sep=""), width=950, height=400)
  autoplot(prcomp(df), data = cdata, colour = 'Cohort', main=base)
  dev.off()
  png(paste("PCA_",base,"_Centre_original.png",sep=""), width=950, height=400)
  autoplot(prcomp(df), data = cdata, colour = 'Centre', main=base)
  dev.off()
  

  write.table(ddata,output_file_name,sep="\t")
  
  ###################################  
  # QN
  ################################### 
  
  QN_normalize(output_file_name,qn_output_file_name)
  
  ###################################  
  # QN PCA
  ################################### 
  
  met <-read.table(qn_output_file_name,header=T,sep="\t")
  cov_file_name <- "/home/natalja/Documents/EMIF/KCL/Covariates_all_vertical_text.txt"
  cov <- read.table(cov_file_name,header = T,sep="\t", stringsAsFactors = FALSE)
  rownames(cov) <- cov$Sample
  cov<-cov[colnames(met),]
  
  e <- met[,names(met) %in% rownames(cov)]
  cov <- cov[rownames(cov) %in% colnames(e),]
  cov <- cov[colnames(e),]
  
  ddata <- as.matrix(e)
  
  scaleShift=abs(min(t(ddata), na.rm = TRUE))+1 
  df = log(t(ddata)+scaleShift)
  cdata <- as.data.frame(df)
  cdata$Diagnosis <- cov$Diagnosis
  cdata$Cohort <- cov$Cohort
  cdata$Centre <- cov$Centre
  cdata$Gender <- cov$Gender
  png(paste("PCA_",base,"_Diagnosis_QN.png",sep=""), width=950, height=400)
  autoplot(prcomp(df), data = cdata, colour = 'Diagnosis', main=base)
  dev.off()
  png(paste("PCA_",base,"_Cohort_QN.png",sep=""), width=950, height=400)
  autoplot(prcomp(df), data = cdata, colour = 'Cohort', main=base)
  dev.off()
  png(paste("PCA_",base,"_Centre_QN.png",sep=""), width=950, height=400)
  autoplot(prcomp(df), data = cdata, colour = 'Centre', main=base)
  dev.off()
}  

