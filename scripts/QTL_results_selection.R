################################################################################
# Subsets metabolomic matrices based on metabolic QTL results
################################################################################
feature_selection_mqtl_per_assay <- function(base,metabolites=TRUE){
  threshold <- 0.01
  dna_matrices_dir<-"./data/dna_matrices/"
  met_matrices_dir<-"./data/metabolomics_matrices/normalised/"
  cov_matrices_dir<-"./data/covariates/"
  
  res <- read.table(paste("mqtl_",base,".txt",sep=""),header=T,sep="\t")
  res_sig <- res[res$FDR<threshold,]

  if (metabolites){
    write.table(unique(res_sig$gene),paste(base,"_met_sig.txt",sep=""),sep="\t",row.names=F,quote=F)
  }
  else {
    write.table(unique(res_sig$SNP),paste(base,"_SNP_sig.txt",sep=""),sep="\t",row.names=F,quote=F)
  }
  
  # Matrices are to big to be processed in R
  # print(paste("grep -wFf ",base,"_SNP_sig.txt ",dna_matrices_dir,"dna_matrix_",base,".tsv > dna_",base,"_sig.txt",sep=""))
  # print(paste("grep -wFf ",base,"_met_sig.txt ",met_matrices_dir,base,".tsv > met_",base,"_sig.txt",sep=""))
  # print(paste("head -1 ",met_matrices_dir,base,".tsv > ",base,"_header.txt",sep=""))
  # print(paste("cat ",base,"_header.txt met_",base,"_sig.txt > met_",base,"_sig2.txt",sep=""))
  # print(paste("cat ",base,"_header.txt dna_",base,"_sig.txt > dna_",base,"_sig2.txt",sep=""))
}

################################################################################
# List of metabolites from significant QTL results, all assays together
################################################################################
feature_selection_mqtl_metabolites <- function(){
  # LC-MS
  base <- "UHPOS"
  expr <- read.table(paste("./results/",base,"_met_sig.txt",sep=""),sep="\t",header=T,row.names=1,stringsAsFactors = F)
  expr_urine_all <-expr[0,]
  expr_urine_all$base <= ""
  for(base in c("UHPOS","URPOS","URNEG")){
    expr <- read.table(paste("./results/",base,"_met_sig.txt",sep=""),sep="\t",header=T,row.names=1,stringsAsFactors = F)
    expr$base <- base
    expr_urine_all <- rbind(expr_urine_all,expr)
  }


  # NMR 
  base <-"NMR_Urine"
  expr_nmr <- read.table(paste("./results/",base,"_met_sig.txt",sep=""),sep="\t",header=T,row.names=1,stringsAsFactors = F)
  expr_nmr$base <- "NMR"

  expr_urine_all <- rbind(expr_urine_all,expr_nmr)
  
  write.table(expr_urine_all,"./results/metabolites.txt",sep="\t",row.names=F,quote=F)
}

################################################################################
# List of SNPs from QTL significant results, all assays together
################################################################################
feature_selection_mqtl_snps <- function(){
  threshold <- 0.01
  # LC-MS
  base <- "UHPOS"
  res <- read.table(paste("./results/mqtl_",base,".txt",sep=""),header=T,sep="\t")
  
  res_urine_all <-res[0,]
  res_urine_all$base <= ""
  for(base in c("UHPOS","URPOS","URNEG")){
    res <- read.table(paste("./results/mqtl_",base,".txt",sep=""),header=T,sep="\t")
    res_sig <- res[res$FDR<threshold,]
    res_sig$SNP <- gsub('_2','',res_sig$SNP)
    res_sig$SNP <- gsub('_1','',res_sig$SNP)
    res_sig$base <- base
    res_urine_all <- rbind(res_urine_all,res_sig)
  }
  
  base <-"NMR_Urine"
  res_nmr <- read.table(paste("./results/mqtl_",base,".txt",sep=""),sep="\t",header=T)
  res_nmr$base <- "NMR"
  res_nmr <- res_nmr[res_nmr$FDR<threshold,]
  res_nmr$SNP <- gsub('_2','',res_nmr$SNP)
  res_nmr$SNP <- gsub('_1','',res_nmr$SNP)
  
  snps_urine_all <- unique(c(unique(res_urine_all$SNP),unique(res_nmr$SNP)))
  
  write.table(snps_urine_all,"./results/snps.txt",sep="\t",row.names=F,quote=F)
}

################################################################################
# SNPs from mQTL and GWAS intersection
################################################################################
snps_intersection <- function(gwas_snps_file){
  gwas_snps <- read.table(gwas_snps_file,sep="\t",header = T,stringsAsFactors = F)
  snps_urine_all <- read.table("./results/snps.txt",sep="\t",header = T,stringsAsFactors = F)
  intersect(gwas_snps$GWAS.SNP,snps_urine_all)
}

