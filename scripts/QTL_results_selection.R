################################################################################
# Subsets metabolomic matrices based on metabolic QTL results
################################################################################
threshold <- 0.01
base <- "UHPOS"

#res <- read.table(paste("/hps/nobackup/ma/natalja/data/mQTL_data/NMR_Results/previous_results/mqtl_",base,".txt",sep=""),header=T,sep="\t")
dna_matrices_dir<-"/hps/nobackup/ma/natalja/data/mQTL_data/dna_matrices/"
met_matrices_dir<-"/hps/nobackup/ma/natalja/data/mQTL_data/met_matrices/normalization_final/" #/NMR/
cov_matrices_dir<-"/hps/nobackup/ma/natalja/data/mQTL_data/covariates/"

res <- read.table(paste("mqtl_",base,".txt",sep=""),header=T,sep="\t")
res_sig <- res[res$FDR<threshold,]
#res_sig$SNP <- gsub('_2','',res_sig$SNP)
#res_sig$SNP <- gsub('_1','',res_sig$SNP)

# Preparation
dim(res_sig)
length(unique(res_sig$SNP))
length(unique(res_sig$gene))

write.table(unique(res_sig$SNP),paste(base,"_sig_SNP.txt",sep=""),sep="\t",row.names=F,quote=F)
write.table(unique(res_sig$gene),paste(base,"_sig_met.txt",sep=""),sep="\t",row.names=F,quote=F)

print(paste("grep -wFf ",base,"_sig_SNP.txt ",dna_matrices_dir,"dna_matrix_",base,".tsv > dna_",base,"_sig.txt",sep=""))
print(paste("grep -wFf ",base,"_sig_met.txt ",met_matrices_dir,base,".tsv > met_",base,"_sig.txt",sep=""))
print(paste("head -1 ",met_matrices_dir,base,".tsv > ",base,"_header.txt",sep=""))
print(paste("cat ",base,"_header.txt met_",base,"_sig.txt > met_",base,"_sig2.txt",sep=""))
print(paste("cat ",base,"_header.txt dna_",base,"_sig.txt > dna_",base,"_sig2.txt",sep=""))

################################################################################
# List of metabolites from QTL results
# List of significant mQTL results from all assays
################################################################################
threshold <- 0.01
drugs <- read.table("/hps/nobackup/ma/natalja/data/mQTL_data/drugs.txt",sep="\t",header = T,stringsAsFactors = F)

# LC-MS
base <- "UHPOS"
expr <- read.table(paste("/hps/nobackup/ma/natalja/data/mQTL_data/results/normalization_final/met_",base,"_sig2.txt",sep=""),sep="\t",header=T,row.names=1,stringsAsFactors = F)

expr_urine_all <-expr[0,]
expr_urine_all$base <= ""
for(base in c("UHPOS","URPOS","URNEG")){
  expr <- read.table(paste("/hps/nobackup/ma/natalja/data/mQTL_data/results/normalization_final/met_",base,"_sig2.txt",sep=""),sep="\t",header=T,row.names=1,stringsAsFactors = F)
  expr$base <- base
  expr_urine_all <- rbind(expr_urine_all,expr)
}


# NMR 
base <-"NMR_Nosey_Urine"
expr_nmr <- read.table(paste("/hps/nobackup/ma/natalja/data/mQTL_data/NMR_Results/normalization_final2/met_",base,"_sig2.txt",sep=""),sep="\t",header=T,row.names=1,stringsAsFactors = F)
expr_nmr$base <- "NMR NOESY"

expr_urine_all <- rbind(expr_urine_all,expr_nmr)

write.table(expr_urine_all,"metabolites.txt",sep="\t",row.names=F,quote=F)

################################################################################
# SNPs from mQTL results
################################################################################
gwas_snps <- read.table("/hps/nobackup/ma/natalja/data/mQTL_data/gwas_snps.txt",sep="\t",header = T,stringsAsFactors = F)
extra_snps <- read.table("/hps/nobackup/ma/natalja/data/mQTL_data/gwas_snps2.txt",sep="\t",header = T,stringsAsFactors = F)
our_gwas <- c(
  "rs429358",
  "rs769449",
  "rs56131196",
  "rs4420638",
  "rs10414043",
  "rs7256200",
  "rs12721051",
  "rs73052335",
  "rs12721046"
)

all_snps <- read.table("/hps/nobackup/ma/natalja/data/mQTL_data/snps.txt",sep="\t",header = T,stringsAsFactors = F)
all_snps$Sample <- gsub('_2','',all_snps$Sample)
all_snps$Sample <- gsub('_1','',all_snps$Sample)
all_snps$Sample <- gsub(':D','',all_snps$Sample)
all_snps$Sample <- gsub(':I','',all_snps$Sample)

threshold <- 0.01
# LC-MS
base <- "UHPOS"
res <- read.table(paste("/hps/nobackup/ma/natalja/data/mQTL_data/results/normalization_final/mqtl_",base,".txt",sep=""),header=T,sep="\t")

res_urine_all <-res[0,]
res_urine_all$base <= ""
for(base in c("UHPOS","URPOS","URNEG")){
  res <- read.table(paste("/hps/nobackup/ma/natalja/data/mQTL_data/results/normalization_final/mqtl_",base,".txt",sep=""),header=T,sep="\t")
  res_sig <- res[res$FDR<threshold,]
  res_sig$SNP <- gsub('_2','',res_sig$SNP)
  res_sig$SNP <- gsub('_1','',res_sig$SNP)
  res_sig$base <- base
  res_urine_all <- rbind(res_urine_all,res_sig)
}

base <-"NMR_Nosey_Urine"
res_nmr <- read.table(paste("/hps/nobackup/ma/natalja/data/mQTL_data/NMR_Results/normalization_final2/mqtl_",base,".txt",sep=""),sep="\t",header=T)
res_nmr$base <- "NMR NOESY"
res_nmr <- res_nmr[res_nmr$FDR<threshold,]
res_nmr$SNP <- gsub('_2','',res_nmr$SNP)
res_nmr$SNP <- gsub('_1','',res_nmr$SNP)

snps_urine_all <- unique(c(unique(res_urine_all$SNP),unique(res_nmr$SNP)))
intersect(gwas_snps$GWAS.SNP,snps_urine_all)

