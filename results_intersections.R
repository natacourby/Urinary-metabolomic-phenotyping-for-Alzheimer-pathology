################################################################################
# Previous studies 
################################################################################
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3446258/
ps <- read.table("/hps/nobackup/ma/natalja/data/mQTL_data/previous_studies.txt",sep="\t",header = T,stringsAsFactors = F)

for(base in c("UHPOS","URPOS","URNEG")){
  main_dir <- "/hps/nobackup/ma/natalja/data/mQTL_data/results/"
  
  output_file_name <- paste(main_dir,"result_",base,"_annot_full.txt",sep="") 
  res_all <- read.table(output_file_name,header=T,sep="\t")
  print(base)
  print(intersect(res_all$SNP, ps$SNP.ID))
  #write.table(res_all[res_all$SNP %in% intersect(res_all$SNP, ps$SNP.ID),],paste(main_dir,"previous_studies_",base,".txt",sep=""),sep="\t")
}

for(base in c("SHPOS","SLPOS","SLNEG")){
  main_dir <- "/hps/nobackup/ma/natalja/data/mQTL_data/results/normalization_final/"
  output_file_name <- paste(main_dir,"result_",base,"_annot_full.txt",sep="") 
  res_all <- read.table(output_file_name,header=T,sep="\t",stringsAsFactors = F)
  print(base)
  print(intersect(res_all$SNP, ps$SNP.ID))
  #write.table(res_all[res_all$SNP %in% intersect(res_all$SNP, ps$SNP.ID),],paste(main_dir,"previous_studies_",base,".txt",sep=""),sep="\t")

}
  
for (base in c("NMR_CPMG_Serum","NMR_Lipidomics_Negative_Serum","NMR_Lipidomics_Positive_Serum","NMR_Lipidomics_Negative_Urine","NMR_Lipidomics_Positive_Urine","NMR_Nosey_Urine"))
{
  main_dir <- "/hps/nobackup/ma/natalja/data/mQTL_data/NMR_Results/normalization_final/"
  output_file_name <- paste(main_dir,"result_",base,"_annot_full.txt",sep="")   
  res_all <- read.table(output_file_name,header=T,sep="\t",stringsAsFactors = F)
  print(base)
  print(intersect(res_all$SNP, ps$SNP.ID))
  
} 


################################################################################
# LC-MS and NMR intersection  
################################################################################

res_serum <-res_all[0,]
for(base in c("SHPOS","SLPOS","SLNEG")){
  main_dir <- "/hps/nobackup/ma/natalja/data/mQTL_data/results/EigenMS_QN_removed_cohort_centre/"
  output_file_name <- paste(main_dir,"result_",base,"_annot_full.txt",sep="")   
  res_all <- read.table(output_file_name,header=T,sep="\t",stringsAsFactors = F)
  res_serum <- rbind(res_serum,res_all)
}

res_serum_nmr <-res_all[0,]
for (base in c("NMR_CPMG_Serum","NMR_Lipidomics_Negative_Serum","NMR_Lipidomics_Positive_Serum"))
{
  main_dir <- "/hps/nobackup/ma/natalja/data/mQTL_data/NMR_Results/normalization_final/"
  output_file_name <- paste(main_dir,"result_",base,"_annot_full.txt",sep="") 
  res_all <- read.table(output_file_name,header=T,sep="\t",stringsAsFactors = F)
  res_serum_nmr <- rbind(res_serum_nmr,res_all)
} 

list_snps <- intersect(unique(res_serum$SNP),unique(res_serum_nmr$SNP))

unique(res_serum[res_serum$SNP %in% list_snps,]$gene)
unique(res_serum_nmr[res_serum_nmr$SNP %in% list_snps,]$gene)
# 94
# 68
res_urine <-res_all[0,]
for(base in c("UHPOS","URPOS","URNEG")){
  main_dir <- "/hps/nobackup/ma/natalja/data/mQTL_data/results/"
  output_file_name <- paste(main_dir,"result_",base,"_annot_full.txt",sep="")   
  res_all <- read.table(output_file_name,header=T,sep="\t",stringsAsFactors = F)
  res_urine <- rbind(res_urine,res_all)
}


res_urine_nmr <-res_all[0,]
for (base in c("NMR_Lipidomics_Negative_Urine","NMR_Lipidomics_Positive_Urine","NMR_Nosey_Urine"))
{
  main_dir <- "/hps/nobackup/ma/natalja/data/mQTL_data/NMR_Results/normalization_final/"
  output_file_name <- paste(main_dir,"result_",base,"_annot_full.txt",sep="") 
  res_all <- read.table(output_file_name,header=T,sep="\t",stringsAsFactors = F)
  res_urine_nmr <- rbind(res_urine_nmr,res_all)
} 
list_snps <- intersect(unique(res_urine$SNP),unique(res_urine_nmr$SNP))

unique(res_urine[res_urine$SNP %in% list_snps,]$gene)
unique(res_urine_nmr[res_urine_nmr$SNP %in% list_snps,]$gene)
# 379
# 261

