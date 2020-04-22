################################################################################
# Previous studies
################################################################################
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3446258/
ps <- read.table("/hps/nobackup/ma/natalja/data/mQTL_data/previous_studies_all.txt",sep="\t",header = T,stringsAsFactors = F)
drugs <- read.table("/hps/nobackup/ma/natalja/data/mQTL_data/drugs.txt",sep="\t",header = T,stringsAsFactors = F)
threshold <- 0.01

#########################################
# Load metabolomic selected matrices 
#########################################
# LC-MS
base <- "UHPOS"
main_dir <- "/hps/nobackup/ma/natalja/data/mQTL_data/results/normalization_final/"
output_file_name <- paste(main_dir,"result_",base,"_annot_full.txt",sep="")
res_all <- read.table(output_file_name,header=T,sep="\t")
expr <- read.table(paste("/hps/nobackup/ma/natalja/data/mQTL_data/results/normalization_final/met_",base,"_sig2.txt",sep=""),sep="\t",header=T,row.names=1,stringsAsFactors = F)

res_urine_all <-res_all[0,]
res_urine_all$base <= ""
expr_urine_all <-expr[0,]
#expr_urine_all$base <= ""
for(base in c("UHPOS","URPOS","URNEG")){
  main_dir <- "/hps/nobackup/ma/natalja/data/mQTL_data/results/normalization_final/"
  output_file_name <- paste(main_dir,"result_",base,"_annot_full.txt",sep="")
  res_all <- read.table(output_file_name,header=T,sep="\t",stringsAsFactors = F)
  res_all$base <- base
  expr <- read.table(paste("/hps/nobackup/ma/natalja/data/mQTL_data/results/normalization_final/met_",base,"_sig2.txt",sep=""),sep="\t",header=T,row.names=1,stringsAsFactors = F)
  #expr$base <- base
  res_urine_all <- rbind(res_urine_all,res_all)
  expr_urine_all <- rbind(expr_urine_all,expr)
}
res_sig_lcms <- res_urine_all[res_urine_all$FDR<=threshold,]
res_sig_lcms$SNP <- gsub('_2','',res_sig_lcms$SNP)
res_sig_lcms$SNP <- gsub('_1','',res_sig_lcms$SNP)

# NMR 
base <-"NMR_Nosey_Urine"
expr_nmr <- read.table(paste("/hps/nobackup/ma/natalja/data/mQTL_data/NMR_Results/normalization_final2/met_",base,"_sig2.txt",sep=""),sep="\t",header=T,row.names=1,stringsAsFactors = F)
res <- read.table(paste("/hps/nobackup/ma/natalja/data/mQTL_data/NMR_Results/normalization_final2/mqtl_",base,".txt",sep=""),header=T,sep="\t")
res$base <- base
res_sig_nmr <- res[res$FDR<=threshold,]
res_sig_nmr$SNP <- gsub('_2','',res_sig_nmr$SNP)
res_sig_nmr$SNP <- gsub('_1','',res_sig_nmr$SNP)

# Remove drugs
expr_lcms <- expr_urine_all[!rownames(expr_urine_all) %in% drugs$Metabolite,]
dim(expr_lcms)
# There are not know drugs in NMR
dim(expr_nmr)

ps$LC.MS <- ps$SNP %in% res_sig_lcms$SNP
ps$NMR <- ps$SNP %in% res_sig_nmr$SNP

write.table(ps,"/hps/nobackup/ma/natalja/data/mQTL_data/previous_results_all.txt",sep="\t",row.names=FALSE)

list_snps <- intersect(unique(res_sig_lcms$SNP),unique(res_sig_nmr$SNP))

#########################################
# Venn Diagrams
#########################################
library(VennDiagram)
library(qqman)
png(paste("VennDiagram_Urine.png",sep=""), width=950, height=400)
grid.newpage()
draw.pairwise.venn(length(unique(res_sig_lcms$SNP)), length(unique(res_sig_nmr$SNP)), length(list_snps), category = c("LC-MS SNPs", "NMR SNPs"),
                   lty = rep("blank", 2), fill = c("light blue", "pink"), alpha = rep(0.5, 2), cat.pos = c(0,0), cat.dist = rep(0.025, 2))
dev.off()

png(paste("VennDiagram_ps_SNPs.png",sep=""), width=950, height=400)
grid.newpage()
draw.triple.venn(length(unique(res_sig_lcms$SNP)), length(unique(res_sig_nmr$SNP)),  length(unique(ps$SNP)),
                 length(intersect(res_sig_lcms$SNP,res_sig_nmr$SNP)), 
                  length(intersect(res_sig_nmr$SNP,ps$SNP)), 
                 length(intersect(res_sig_lcms$SNP,ps$SNP)),
                 length(intersect(intersect(res_sig_lcms$SNP,ps$SNP),res_sig_nmr$SNP)),
                 category = c("LC-MS SNPs","NMR SNPs", "Previous studies"),
                 fill = c("light blue", "pink", "green"),
                 lty = "blank",
                 cex = 2,
                 cat.cex = 2)

#draw.pairwise.venn( length(unique(ps$SNP)), length(unique(res_sig_nmr$SNP)), length(list_snps), category = c("LC-MS SNPs", "NMR SNPs"),
 #                  lty = rep("blank", 2), fill = c("light blue", "pink"), alpha = rep(0.5, 2), cat.pos = c(0,0), cat.dist = rep(0.025, 2))
dev.off()

png(paste("VennDiagram_ps_genes.png",sep=""), width=950, height=400)
grid.newpage()
draw.triple.venn(489,38,65,
                 27, 
                 10, 
                 21,
                 10,
                 category = c("LC-MS genes","NMR genes", "Previous studies"),
                 fill = c("light blue", "pink", "green"),
                 lty = "blank",
                 cex = 2,
                 cat.cex = 2)

#draw.pairwise.venn( length(unique(ps$SNP)), length(unique(res_sig_nmr$SNP)), length(list_snps), category = c("LC-MS SNPs", "NMR SNPs"),
#                  lty = rep("blank", 2), fill = c("light blue", "pink"), alpha = rep(0.5, 2), cat.pos = c(0,0), cat.dist = rep(0.025, 2))
dev.off()
# Number of tested SNPs in NMR and LC-MS
length(unique(res_urine_all$SNP)) 
# 1358611
length(unique(res$SNP)) 
# 694279
length(intersect(res$SNP,res_urine_all$SNP)) 
# 186398


#########################################
# Possible names for metabolites
#########################################
base <- "UHPOS"
main_dir <- "/hps/nobackup/ma/natalja/data/mQTL_data/results/normalization_final/"
output_file_name <- paste(main_dir,"result_",base,"_annot_full.txt",sep="")
res_all <- read.table(output_file_name,header=T,sep="\t")

# LC-MS Urine
res_all$subset <- ""
res_urine <- res_all[0,]
for(base in c("UHPOS","URPOS","URNEG")){
  main_dir <- "/hps/nobackup/ma/natalja/data/mQTL_data/results/normalization_final/"
  output_file_name <- paste(main_dir,"result_",base,"_annot_full.txt",sep="")
  res_all <- read.table(output_file_name,header=T,sep="\t",stringsAsFactors = FALSE)
  res_all$subset <- base
  print(base)
  print(intersect(res_all$SNP, ps$SNP.ID))
  res_urine <- rbind(res_urine,res_all[res_all$SNP %in% intersect(res_all$SNP,ps$SNP.ID),])
}
res_urine <- merge(res_urine,ps,by.x="SNP",by.y="SNP.ID",all.x=T,all.y=F)
res <- data.frame(SNP=res_urine$SNP, metabolite=res_urine$gene, subset = res_urine$subset, FDR = res_urine$FDR, possible_name=res_urine$Metabolite, biofluid = res_urine$Biofluid, publication=res_urine$Publication)
res <- unique(res)
write.table(res,"previous_results_ms.txt",sep="\t",row.names=FALSE)

# NMR Urine
res_all$subset <- ""
res_urine_nmr <- res_all[0,]
for (base in c("NMR_Nosey_Urine"))
{
  main_dir <- "/hps/nobackup/ma/natalja/data/mQTL_data/NMR_Results/normalization_final2/"
  output_file_name <- paste(main_dir,"result_",base,"_annot_full.txt",sep="")
  res_all <- read.table(output_file_name,header=T,sep="\t",stringsAsFactors = F)
  res_all$subset <- base
  print(base)
  print(intersect(res_all$SNP, ps$SNP.ID))
  res_urine_nmr <- rbind(res_urine_nmr,res_all[res_all$SNP %in% intersect(res_all$SNP,ps$SNP.ID),])
}
res_urine_nmr <- merge(res_urine_nmr,ps,by.x="SNP",by.y="SNP.ID",all.x=T,all.y=F)
res <- data.frame(SNP=res_urine_nmr$SNP, metabolite=res_urine_nmr$gene, subset = res_urine_nmr$subset, FDR = res_urine_nmr$FDR, possible_name=res_urine_nmr$Metabolite, biofluid = res_urine_nmr$Biofluid)
res <- unique(res)
write.table(res,"previous_results_nmr.txt",sep="\t",row.names=FALSE)

#########################################
# LC-MS and NMR results intersection
#########################################
library(VennDiagram)
library(qqman)
main_dir <- "/hps/nobackup/ma/natalja/data/mQTL_data/results/"
output_file_name <- paste(main_dir,"result_SHPOS_annot_full.txt",sep="")
res_all <- read.table(output_file_name,header=T,sep="\t",stringsAsFactors = F)

# Serum
##########
res_serum <-res_all[0,]
for(base in c("SHPOS","SLPOS","SLNEG")){
  main_dir <- "/hps/nobackup/ma/natalja/data/mQTL_data/results/normalization_final/"
  output_file_name <- paste(main_dir,"result_",base,"_annot_full.txt",sep="")
  res_all <- read.table(output_file_name,header=T,sep="\t",stringsAsFactors = F)
  res_serum <- rbind(res_serum,res_all)
}

res_serum_nmr <-res_all[0,]
for (base in c("NMR_CPMG_Serum","NMR_Nosey_Serum"))
{
  main_dir <- "/hps/nobackup/ma/natalja/data/mQTL_data/NMR_Results/normalization_final2/"
  output_file_name <- paste(main_dir,"result_",base,"_annot_full.txt",sep="")
  res_all <- read.table(output_file_name,header=T,sep="\t",stringsAsFactors = F)
  res_serum_nmr <- rbind(res_serum_nmr,res_all)
}

# List of SNPs that are the same between LC-MS and NMR significant mQTL associations
list_snps <- intersect(unique(res_serum$SNP),unique(res_serum_nmr$SNP))
write.table(list_snps,"serum_common_snps.txt",sep="\t",row.names=FALSE)

print(length(unique(res_serum$SNP)))
print(length(unique(res_serum_nmr$SNP)))
print(length(list_snps))

#########################################
# Venn Diagram of SNPs between LC-MS and NMR results
#########################################
png(paste("VennDiagram_Serum.png",sep=""), width=950, height=400)
grid.newpage()
draw.pairwise.venn(length(unique(res_serum$SNP)), length(unique(res_serum_nmr$SNP)), length(list_snps), category = c("LC-MS SNPs", "NMR SNPs"),
                   lty = rep("blank", 2), fill = c("light blue", "pink"), alpha = rep(0.5, 2), cat.pos = c(0,0), cat.dist = rep(0.025, 2))
dev.off()


#########################################
# Manhattan plot for SNPs in intersection
#########################################
res_serum_all <-res_all[0,]
for (base in c("NMR_Nosey_Serum","NMR_CPMG_Serum"))
{
  main_dir <- "/hps/nobackup/ma/natalja/data/mQTL_data/NMR_Results/normalization_final2/"
  output_file_name <- paste(main_dir,"mqtl_annot_",base,".txt",sep="")
  res_all <- read.table(output_file_name,header=T,sep="\t",stringsAsFactors = F)
  res_serum_all <- rbind(res_serum_all,res_all)
}

# FDR - qqman manhattan plot 
res_plot<-data.frame(SNP=res_serum_all$SNP,CHR=as.character(res_serum_all$chr_name),P=res_serum_all$FDR,
                     BP=res_serum_all$chrom_start,stringsAsFactors=FALSE)
res_plot$CHR[res_plot$CHR == "X"] <- "23"
res_plot$CHR <- as.numeric(res_plot$CHR)

png(paste("serum_SNPs_in_intersection_manhattan_FDR_mqtl_",base,"_nmr.png",sep=""), width=950, height=400)
manhattan(res_plot, ylim=c(0,35), main = "Urine NMR",
          suggestiveline=-log10(1e-2),
          genomewideline=FALSE,
          highlight=intersect(res_plot$SNP,list_snps))
dev.off()

res_serum_all <-res_all[0,]
for(base in c("SHPOS","SLPOS","SLNEG")){
  main_dir <- "/hps/nobackup/ma/natalja/data/mQTL_data/results/"
  output_file_name <- paste(main_dir,"mqtl_annot_",base,".txt",sep="")
  res_all <- read.table(output_file_name,header=T,sep="\t",stringsAsFactors = F)
  res_serum_all <- rbind(res_serum_all,res_all)
}

# FDR - qqman manhattan plot 
res_plot<-data.frame(SNP=res_serum_all$SNP,CHR=as.character(res_serum_all$chr_name),P=res_serum_all$FDR,
                     BP=res_serum_all$chrom_start,stringsAsFactors=FALSE)
res_plot$CHR[res_plot$CHR == "X"] <- "23"
res_plot$CHR <- as.numeric(res_plot$CHR)

png(paste("serum_SNPs_in_intersection_manhattan_FDR_mqtl_",base,"_ms.png",sep=""), width=950, height=400)
manhattan(res_plot, ylim=c(0,35), main = "Urine LC-MS",
          suggestiveline=-log10(1e-2),
          genomewideline=FALSE,
          highlight=intersect(res_plot$SNP,list_snps))
dev.off()

#########################################
# LC-MS, NMR metabolites paired
#########################################
# Serum
met_ms <- unique(res_serum[res_serum$SNP %in% list_snps,]$gene)
met_nmr <- unique(res_serum_nmr[res_serum_nmr$SNP %in% list_snps,]$gene)

res_ms <- res_serum[res_serum$gene %in% met_mc,]

res_nmr <- res_serum_nmr[res_serum_nmr$gene %in% met_nmr,]

res <- merge(res_ms,res_nmr, by="SNP", all=FALSE)
write.table(res,"serum_ms_nmr_intersection.txt",sep="\t")
# unique pairs mc nmr metabolites from Excel
pairs <- read.table("serum_ms_pairs.txt",sep="\t",header=T, stringsAsFactors=F)
p2<-pairs[duplicated(pairs$met_ms),]
write.table(p2[order(p2$met_mc),],"../serum_nmr_pairs.txt",sep="\t")

# Urine
res_urine <-res_all[0,]
for(base in c("UHPOS","URPOS","URNEG")){
  main_dir <- "/hps/nobackup/ma/natalja/data/mQTL_data/results/normalization_final/"
  output_file_name <- paste(main_dir,"result_",base,"_annot_full.txt",sep="")
  res_all <- read.table(output_file_name,header=T,sep="\t",stringsAsFactors = F)
  res_urine <- rbind(res_urine,res_all)
}

res_urine_nmr <-res_all[0,]
for (base in c("NMR_Nosey_Urine"))
{
  main_dir <- "/hps/nobackup/ma/natalja/data/mQTL_data/NMR_Results/normalization_final2/"
  output_file_name <- paste(main_dir,"result_",base,"_annot_full.txt",sep="")
  res_all <- read.table(output_file_name,header=T,sep="\t",stringsAsFactors = F)
  res_urine_nmr <- rbind(res_urine_nmr,res_all)
}
list_snps <- intersect(unique(res_urine$SNP),unique(res_urine_nmr$SNP))
write.table(list_snps,"urine_common_snps.txt",sep="\t",row.names=FALSE)

print(length(unique(res_urine$SNP)))
print(length(unique(res_urine_nmr$SNP)))
print(length(list_snps))

met_ms <- unique(res_urine[res_urine$SNP %in% list_snps,]$gene)
met_nmr <- unique(res_urine_nmr[res_urine_nmr$SNP %in% list_snps,]$gene)

res_ms <- res_urine[res_urine$gene %in% met_mc,]

res_nmr <- res_urine_nmr[res_urine_nmr$gene %in% met_nmr,]

res_ms_subset<-data.frame(SNP=res_mc$SNP,met_mc=res_mc$gene,chr=res_mc$chr_name,position=res_mc$chrom_start,gene=res_mc$hgnc_symbol,ensembl_id=res_mc$ensembl_gene_st)
res_nmr_subset<-data.frame(SNP=res_nmr$SNP,met_nmr=res_nmr$gene)

res <- merge(res_ms_subset,res_nmr_subset, by="SNP", all=FALSE)
write.table(res,"urine_ms_nmr_intersection.txt",sep="\t")

write.table(res_urine[res_urine$SNP %in% intersect(res_urine$SNP, ps$SNP.ID),],"urine_ms_ps.txt",sep="\t")
write.table(res_urine_nmr[res_urine_nmr$SNP %in% intersect(res_urine_nmr$SNP, ps$SNP.ID),],"urine_nmr_ps.txt",sep="\t")


pairs <- data.frame(met_ms=res$met_ms,met_nmr=res$met_nmr)
pairs <- unique(pairs)

p2<-pairs[duplicated(pairs$met_ms),]
write.table(p2[order(p2$met_mc),],"urine_ms_pairs.txt",sep="\t",row.names=FALSE)

p2<-pairs[duplicated(pairs$met_nmr),]
write.table(p2[order(p2$met_nmr),],"urine_nmr_pairs.txt",sep="\t",row.names=FALSE)