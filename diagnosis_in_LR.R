################################################################################
# Linear Regression with diagnosis in a model 
################################################################################

threshold <- 0.01

res <- read.table(paste("/hps/nobackup/ma/natalja/data/mQTL_data/NMR_Results/mqtl_",base,".txt",sep=""),header=T,sep="\t")
res <- read.table(paste("mqtl_",base,".txt",sep=""),header=T,sep="\t")
res_sig <- res[res$FDR<threshold,]
res_sig$SNP <- gsub('_2','',res_sig$SNP)
res_sig$SNP <- gsub('_1','',res_sig$SNP)

dim(res_sig)
length(unique(res_sig$SNP))
length(unique(res_sig$gene))

write.table(unique(res_sig$SNP),paste(base,"_sig_SNP.txt",sep=""),sep="\t",row.names=F,quote=F)
write.table(unique(res_sig$gene),paste(base,"_sig_met.txt",sep=""),sep="\t",row.names=F,quote=F)

# grep -wFf sig_SNP.txt ../dna_matrices/dna_matrix_NMR_Lipidomics_Positive_Urine.tsv > dna_sig.txt
# grep -wFf sig_met.txt ../met_matrices/NMR/NMR_Lipidomics_Positive_Urine.tsv > met_sig.txt
# head -1  ../met_matrices/NMR/NMR_Lipidomics_Positive_Urine.tsv > header.txt
# cat header.txt met_sig.txt > met_sig2.txt 
# cat header.txt dna_sig.txt > dna_sig2.txt
expr <- read.table("met_sig2.txt",sep="\t",header=T,row.names=1,stringsAsFactors = F)
snps <- read.table("dna_sig2.txt",sep="\t",header=T,row.names=1,stringsAsFactors = F)
cvrt = read.table("/hps/nobackup/ma/natalja/data/mQTL_data/covariates/Covariates_NMR_Lipidomics_Positive_Urine.txt",sep="\t",header=T,row.names=1)



cov <-  read.table("/hps/nobackup/ma/natalja/data/mQTL_data/pheno.txt",sep="\t",header=T)
met_qn <-read.table(paste("/hps/nobackup/ma/natalja/data/mQTL_data/met_matrices/NMR/",base,".tsv",sep=""),
                    header=T,sep="\t")
res <- read.table(paste("/hps/nobackup/ma/natalja/data/mQTL_data/NMR_Results/result_",base,"_annot_full.txt",sep=""),
                  header=T,sep="\t",stringsAsFactors = F)

res <- read.table(paste("result_",base,"_annot_full.txt",sep=""),
                  header=T,sep="\t",stringsAsFactors = F)

rownames(cov) <- cov$FID
cov <-cov[!is.na(cov$Diagnosis),]
cov <- cov[cov$Diagnosis %in% c("ADC","CTL"),]

rownames(met_qn) <- met_qn$id
met_qn <- met_qn[,-1]

e <- met_qn[rownames(met_qn) %in% res$gene,names(met_qn) %in% rownames(cov)]
cov <- cov[rownames(cov) %in% colnames(e),]
cov <- cov[colnames(e),]

e <- t(e)
e <- as.data.frame(e)

w_t <- function(x){
  df <- structure(list(Diagnosis = cov$Diagnosis, 
                       Value = x), .Names = c("Diagnosis","Value"), class = "data.frame"); 
  t <- wilcox.test(formula = Value ~ Diagnosis, data = df);
  return(t$p.value)}

res <- apply(e, 2, w_t)
res_FDR <- p.adjust(res,"fdr")
e_imp <- e[,res_FDR<0.01] # important metabolites 

write.table(colnames(e_imp),paste("/hps/nobackup/ma/natalja/data/mQTL_data/NMR_Results/diagnosis_different_metabolites_",base,".txt",sep=""),sep="\t")

e$diagnosis <- cov$Diagnosis
write.table(e,paste("/home/natalja/Documents/EMIF/KCL/mQTL/mQTL_data/results/metabolites_",base,".txt",sep=""),sep="\t")
