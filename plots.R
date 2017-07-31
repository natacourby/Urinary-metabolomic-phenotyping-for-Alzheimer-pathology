################################################################################
# Histogram and manhattan plots for all p-values
################################################################################
args <- commandArgs(trailingOnly=TRUE)
base <- args[1]

source("/hps/nobackup/ma/natalja/data/mQTL_data/scripts/mqtl.R")
library(biomaRt)
library(qqman)
snp.db <- useMart("ENSEMBL_MART_SNP", dataset="hsapiens_snp")

output_file_name <- paste("mqtl_",base,".txt",sep="") 
res <- read.table(output_file_name,header=T,sep="\t")

res$SNP <- gsub('_2','',res$SNP)
res$SNP <- gsub('_1','',res$SNP)

png(paste("histogram_p_values_",base,".png",sep=""))
hist(res$p.value[!is.na(res$p.value)],breaks = 50,xlab = "P-values",col = "paleturquoise3",main = paste(base,sep=""))
dev.off()

res_split <- split(res, (seq(nrow(res))-1) %/% 20000) 
res_all <- res[0,]
for(i in 1:length(res_split)){
  res_part <- res_split[[i]]
  nt.biomart <- getBM(c("refsnp_id","chr_name","chrom_start"),
                      filters="snp_filter",
                      values=res_part$SNP,
                      mart=snp.db)
  res_all_annot <- merge(res_part,nt.biomart,by.x="SNP",by.y="refsnp_id",all.x=T,all.y=F)
  res_all_annot <- res_all_annot[-grep("CHR_", res_all_annot$chr_name),]
  res_all<- rbind(res_all,res_all_annot)
}

res_all <- res_all[!is.na(res_all$chr_name) & res_all$chr_name!="",]
res_all <- res_all[!is.na(res_all$p.value) & res_all$p.value!="",]
res_all <- res_all[!is.na(res_all$chrom_start) & res_all$chrom_start!="",]
res_all <- res_all[!is.na(res_all$SNP) & res_all$SNP!="",]
write.table(res_all,paste("mqtl_annot_",base,".txt",sep="")
            ,sep="\t",row.names = FALSE)

# FDR - qqman manhattan plot 
res_plot<-data.frame(SNP=res_all$SNP,CHR=as.character(res_all$chr_name),P=res_all$FDR,
                     BP=res_all$chrom_start,stringsAsFactors=FALSE)
res_plot$CHR[res_plot$CHR == "X"] <- "23"
res_plot$CHR <- as.numeric(res_plot$CHR)

png(paste("manhattan_FDR_mqtl_",base,".png",sep=""), width=950, height=400)
manhattan(res_plot, ylim=c(0,35), main = paste(base,sep=""),suggestiveline=-log10(1e-2),genomewideline=FALSE)
dev.off()

# raw p-values manhattan.plot 
res_plot<-data.frame(SNP=res_all$SNP,CHR=as.character(res_all$chr_name),P=res_all$p.value,
                     BP=res_all$chrom_start,stringsAsFactors=FALSE)
res_plot$CHR[res_plot$CHR == "X"] <- "23"
res_plot$CHR <- as.numeric(res_plot$CHR)

png(paste("manhattan_pvalues_mqtl_",base,".png",sep=""), width=950, height=400)
manhattan.plot(res_plot$CHR,res_plot$BP,res_plot$P, ylim=c(0,35), main = paste(base,sep=""))
dev.off()





