##################################################################################################
# Feature selection methods
# - Datasets preparation for classification
# - Linear Regression with diagnosis in a model
##################################################################################################


# Datasets preparation for classification
feature_selection_data_preparation <- function(set, classifier){
  expr1 <- read.table("./results/UHPOS_met_sig.txt",sep="\t",header=T,row.names=1,stringsAsFactors = F)
  expr2 <- read.table("./results/URPOS_met_sig.txt",sep="\t",header=T,row.names=1,stringsAsFactors = F)
  expr3 <- read.table("./results/URNEG_met_sig.txt",sep="\t",header=T,row.names=1,stringsAsFactors = F)
  expr_nmr <- read.table("./results/NMR_Urine_met_sig.txt",sep="\t",header=T,row.names=1,stringsAsFactors = F)
  
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
  
  # A: Metabolites only and D: Metabolites and covariates
  if (set=="A" | set=="D"){
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
  }
  # A: Metabolites only
  if (set=="A"){
    # Classifier I: AD/CTL/cMCI/sMCI
    if (classifier=="I"){
      drops <- c("Gender","Centre","Age","Diagnosis2","Diagnosis3")
      df_result <- df[,!(names(df) %in% drops)]
      names(df_result)[names(df_result) == 'Diagnosis1'] <- 'Diagnosis'
    }
    
    # Classifier II: AD+cMCI/CTL+sMCI
    if (classifier=="II"){
      drops <- c("Gender","Centre","Age","Diagnosis1","Diagnosis3")
      df_result <- df[,!(names(df) %in% drops)]
      names(df_result)[names(df_result) == 'Diagnosis2'] <- 'Diagnosis'
    }
    
    # Classifier III: AD/CTL
    if (classifier=="III"){
      drops <- c("Gender","Centre","Age","Diagnosis2","Diagnosis3")
      df_result <- df[cvrt$Diagnosis %in% c("CTL","ADC"),!(names(df) %in% drops)]
      names(df_result)[names(df_result) == 'Diagnosis1'] <- 'Diagnosis'
    }
  }
  
  # D: Metabolites and covariates
  if (set=="D"){
    # Classifier I: AD/CTL/cMCI/sMCI
    if (classifier=="I"){
      drops <- c("Diagnosis2","Diagnosis3")
      df_result <- df[,!(names(df) %in% drops)]
      names(df_result)[names(df_result) == 'Diagnosis1'] <- 'Diagnosis'
    }
    
    # Classifier II: AD+cMCI/CTL+sMCI
    if (classifier=="II"){
      drops <- c("Diagnosis1","Diagnosis3")
      df_result <- df[,!(names(df) %in% drops)]
      names(df_result)[names(df_result) == 'Diagnosis2'] <- 'Diagnosis'
    }
    
    # Classifier III: AD/CTL
    if (classifier=="III"){
      drops <- c("Diagnosis2","Diagnosis3")
      df_result <- df[cvrt$Diagnosis %in% c("CTL","ADC"),!(names(df) %in% drops)]
      names(df_result)[names(df_result) == 'Diagnosis1'] <- 'Diagnosis'
    }
    
    # Classifier IV: cMCI/sMCI
    if (classifier=="IV"){
      drops <- c("Diagnosis2","Diagnosis3")
      df_result <- df[cvrt$Diagnosis %in% c("cMCI","sMCI"),!(names(df) %in% drops)]
      names(df_result)[names(df_result) == 'Diagnosis1'] <- 'Diagnosis'
    }
  }
  
  # B: Metabolites and SNPs, C: Metabolites, SNPs and covariates
  if (set=="B" | set=="C"){
    dna_all <- read.table("./results/snps.txt",sep="\t",header=T,row.names=1,stringsAsFactors = F)
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
    
    if (set=="B"){
      # Classifier I: AD/CTL/cMCI/sMCI
      if (classifier=="I"){
        drops <- c("Gender","Centre","Age","Diagnosis2","Diagnosis3")
        df_result <- df[,!(names(df) %in% drops)]
        names(df_result)[names(df_result) == 'Diagnosis1'] <- 'Diagnosis'
      }
      # Classifier II: AD+cMCI/CTL+sMCI
      if (classifier=="II"){
        drops <- c("Gender","Centre","Age","Diagnosis1","Diagnosis3")
        df_result <- df[,!(names(df) %in% drops)]
        names(df_result)[names(df_result) == 'Diagnosis2'] <- 'Diagnosis'
      }
      # Classifier III: AD/CTL
      if (classifier=="III"){
        drops <- c("Gender","Centre","Age","Diagnosis2","Diagnosis3")
        df_result <- df[cvrt$Diagnosis %in% c("CTL","ADC"),!(names(df) %in% drops)]
        names(df_result)[names(df_result) == 'Diagnosis1'] <- 'Diagnosis'
      }
    }
    
    if (set=="C"){
      # Classifier I: AD/CTL/cMCI/sMCI
      if (classifier=="I"){
        drops <- c("Diagnosis2","Diagnosis3")
        df_result <- df[,!(names(df) %in% drops)]
        names(df_result)[names(df_result) == 'Diagnosis1'] <- 'Diagnosis'
      }
      # Classifier II: AD+cMCI/CTL+sMCI
      if (classifier=="II"){
        drops <- c("Diagnosis1","Diagnosis3")
        df_result <- df[,!(names(df) %in% drops)]
        names(df_result)[names(df_result) == 'Diagnosis2'] <- 'Diagnosis'
      }
      # Classifier III: AD/CTL
      if (classifier=="III"){
        drops <- c("Diagnosis2","Diagnosis3")
        df_result <- df[cvrt$Diagnosis %in% c("CTL","ADC"),!(names(df) %in% drops)]
        names(df_result)[names(df_result) == 'Diagnosis1'] <- 'Diagnosis'
      }
    }
  }
  
  df_result$Diagnosis<- as.factor(df_result$Diagnosis)
  names(df_result) <- make.names(names(df_result))
  return (df_result)
}


# Linear Regression with diagnosis in a model 
feature_selection_LR <- function(base){
  expr <- read.table(paste("./results/",base,"_met_sig.txt",sep=""),sep="\t",header=T,row.names=1,stringsAsFactors = F)
  snps <- read.table(paste("./results/",base,"_SNP_sig.txt",sep=""),sep="\t",header=T,row.names=1,stringsAsFactors = F)
  cvrt = read.table("./data/covariates/Covariates_diagnosis.txt",sep="\t",header=T,row.names=1)
  
  cvrt <- cvrt[,colnames(cvrt) %in% colnames(expr),]
  cvrt <- cvrt[,colnames(expr)]
  
  res_sig$s2 <- 0
  res_sig$Diagnosis_ADC <- 0
  res_sig$Diagnosis_CTL <- 0
  result <-res_sig[0,]
  r <- 1
  
  for(index in 1:nrow(res_sig)){
    i <- res_sig[index,]
    s1 <- as.numeric(snps[as.character(i$SNP),])
    e1 <- as.numeric(expr[as.character(i$gene),])
    
    s2<-s1[!is.na(s1)]
    e2<-e1[!is.na(s1)]
    cvrt2<-cvrt[,!is.na(s1)]
    age <- as.numeric(cvrt2[1,])
    age[is.na(age)] <- 75
    
    gender <- as.numeric(cvrt2[2,])
    CentreLodz <- as.numeric(cvrt2[3,])
    CentrePerugia <- as.numeric(cvrt2[4,])
    CentreThessaloniki <- as.numeric(cvrt2[5,])
    CentreToulouse <- as.numeric(cvrt2[6,])
    Genotype_batch1<- as.numeric(cvrt2[7,])
    Genotype_batch2 <- as.numeric(cvrt2[8,])
    Diagnosis_ADC <- as.numeric(cvrt2[9,])
    Diagnosis_cMCI <- as.numeric(cvrt2[10,])
    Diagnosis_CTL <- as.numeric(cvrt2[11,])
    
    lm2 = lm(e2~s2+age+gender+CentreLodz+CentrePerugia+CentreThessaloniki+CentreToulouse+Genotype_batch1+Genotype_batch2+Diagnosis_ADC+Diagnosis_cMCI+Diagnosis_CTL)
    
    res <- data.frame(summary(lm2)$coef)
    significant_coefficient_names <- paste(rownames(res[res[,4]<=0.01,]), collapse="_")
    if (grepl("Diagnosis_ADC",   significant_coefficient_names) || grepl("Diagnosis_CTL",   significant_coefficient_names)){
      i$s2 <- res[2,4]
      i$Diagnosis_ADC <- res[11,4]
      i$Diagnosis_CTL <- res[13,4]
      result[r,] <- i
      r <- r + 1
      
      i$SNP <- gsub('_2','',i$SNP)
      i$SNP <- gsub('_1','',i$SNP)
      
      lm3 = lm(e2~s2)
      # png(paste("./plots/",base,"_snp_",as.character(i$SNP),"_",gsub('/','',as.character(i$gene)),".png",sep=""), width=950, height=400)
      # plot(e2 ~ jitter(s2),
      #      col=(s2+1),xaxt="n",xlab="Genotype",ylab="Expression")
      # axis(1,at=c(0:2),labels=c("AA","Aa","aa"))
      # lines(lm3$fitted ~ s2,type="b",pch=15,col="darkgrey")
      # dev.off()

    }
  }

  result$SNP <- gsub('_2','',result$SNP)
  result$SNP <- gsub('_1','',result$SNP)
  
  library(biomaRt)
  # Annotation
  snp.db <- useMart("ENSEMBL_MART_SNP", dataset="hsapiens_snp")
  nt.biomart <- getBM(c("refsnp_id","chr_name","ensembl_gene_stable_id"),
                      filters="snp_filter",
                      values=result$SNP,
                      mart=snp.db)
  res_annot <- merge(result,nt.biomart,by.x="SNP",by.y="refsnp_id",all.x=T,all.y=F)
  res_annot <- res_annot[-grep("CHR_", res_annot$chr_name),]
  genes <- unique(nt.biomart$ensembl_gene_stable_id)
  ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
  genes_annot <- getBM(attributes=c('ensembl_gene_id','hgnc_symbol','chromosome_name','description'), 
                       filters = 'ensembl_gene_id', values = genes, mart = ensembl)
  
  
  res_annot <- merge(res_annot,genes_annot,by.x="ensembl_gene_stable_id",by.y="ensembl_gene_id",all.x=T,all.y=F)
  write.table(res_annot,paste("./results/LR_diagnosis_result_annot_",base,".txt",sep=""),sep="\t")
  
  write.table(result,paste("./results/LR_diagnosis_result_",base,".txt",sep=""),sep="\t")
}


