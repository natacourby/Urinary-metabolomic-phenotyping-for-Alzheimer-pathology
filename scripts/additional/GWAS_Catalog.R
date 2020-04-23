################################################################################
# GWAS CATALOG
################################################################################
library(gwascat)

# Load mQTL results
base <- "UHPOS"
res <- read.table(paste("./results/result_",base,"_annot_full.txt",sep=""),header=T,sep="\t")
res_all <- res[0,]
for(base in c("UHPOS","URPOS","URNEG","NMR_Urine")){
  res <- read.table(paste("./results/result_",base,"_annot_full.txt",sep=""),header=T,sep="\t")
  res$SNP <- gsub('_2','',res$SNP)
  res$SNP <- gsub('_1','',res$SNP)
  res$base <- base
  res_all <- rbind(res_all,res)
}
res_genes <- res_all[!is.na(res_all$chr_name) & res_all$chr_name!="",]
res_genes<-res_genes[-grep("CHR_", res_genes$chr_name),]

# Load GWAS Catalog from downloaded file 
# For publication we've used this version: ftp://ftp.ebi.ac.uk/pub/databases/gwas/releases/2018/03/14/gwas-catalog-associations.tsv
gwas <- read.csv("./data/gwas-catalog-associations.tsv",sep="\t")
gwas <- gwas[nchar(as.character(gwas$CHR_ID))<=2,]
gwas$CHR_ID <- as.character(gwas$CHR_ID)
gwas$CHR_POS <- as.character(gwas$CHR_POS)
gwas[nchar(gwas$CHR_POS)==0,]$CHR_POS=1
gwas[nchar(gwas$CHR_ID)==0,]$CHR_ID=0
gwas$CHR_POS <- as.numeric(gwas$CHR_POS)

gr_gwas= GRanges(seqnames=gwas$CHR_ID,IRanges(gwas$CHR_POS, width=1))
mcols(gr_gwas)$'DISEASE/TRAIT'=as.character(gwas$DISEASE.TRAIT)
mcols(gr_gwas)$SNPS=as.character(gwas$SNPS)
mcols(gr_gwas)$'P-VALUE'=gwas$P.VALUE
mcols(gr_gwas)$PVALUE_MLOG=gwas$PVALUE_MLOG

mcols(gr_gwas)$GENE=as.character(gwas$MAPPED_GENE)
genome(gr_gwas)="GRCh38"
ebicat <- new("gwaswloc", gr_gwas,extractDate="2018-03-13")

#########################################
# Select traits related to Azheimer's disease
#########################################
selected_traits <- unique(c(
  unique(getTraits(ebicat))[grep("Alzh",unique(getTraits(ebicat)))],
  unique(getTraits(ebicat))[grep("APOE",unique(getTraits(ebicat)))],
  unique(getTraits(ebicat))[grep("tau",unique(getTraits(ebicat)))],
  unique(getTraits(ebicat))[grep("Tau",unique(getTraits(ebicat)))],
  unique(getTraits(ebicat))[grep("Parkin",unique(getTraits(ebicat)))],
  unique(getTraits(ebicat))[grep("Huntingt",unique(getTraits(ebicat)))],
  unique(getTraits(ebicat))[grep("Amyotrophic lateral sclerosis",unique(getTraits(ebicat)))],
  unique(getTraits(ebicat))[grep("Schizophren",unique(getTraits(ebicat)))],
  unique(getTraits(ebicat))[grep("dementia",unique(getTraits(ebicat)))],
  unique(getTraits(ebicat))[grep("Dementia",unique(getTraits(ebicat)))],
  unique(getTraits(ebicat))[grep("cognitive",unique(getTraits(ebicat)))],
  unique(getTraits(ebicat))[grep("Cognitive",unique(getTraits(ebicat)))]))

ebicat38 <- ebicat

alz_traits <- subsetByTraits(ebicat38,selected_traits)


# Direct mapping by SNPs
intr = ebicat38[ intersect(getRsids(ebicat38), unique(res_genes$SNP)) ]
sort(table(getTraits(intr)), decreasing=TRUE)

subsetByTraits(intr,"Schizophrenia")@elementMetadata$SNPS
subsetByTraits(intr,"Amyotrophic lateral sclerosis")@elementMetadata$SNPS

#met_list <- unique(res_genes[res_genes$SNP %in% subsetByTraits(intr,"Amyotrophic lateral sclerosis")@elementMetadata$SNPS,]$gene)

write.table(sort(table(getTraits(intr)), decreasing=TRUE),"./results/gwas_traits_by_snps.txt",sep="\t")
write.table(intr,"./results/gwas_intersection_details_by_snps.txt")


# Mapping by genomic regions
gr6.0 = GRanges(seqnames=res_genes$chr_name,
                IRanges(res_genes$chrom_start, width=1))
mcols(gr6.0)$rsid = as.character(res_genes$SNP)
ag = function(x) as(x, "GRanges")
ovraw = suppressWarnings(subsetByOverlaps(ag(ebicat38), gr6.0))
length(ovraw)
# average human gene length 8446 bp sd - 7124
ovaug = suppressWarnings(subsetByOverlaps(ag(ebicat38+15570), gr6.0))
length(ovaug)
rawrs = mcols(ovraw)$SNPS
augrs = mcols(ovaug)$SNPS
ebicat38[augrs]
getTraits(ebicat38[augrs])
intr_by_region = ebicat38[augrs]


intr_alz <- subsetByTraits(intr,selected_traits)

intr_alz_by_region <- subsetByTraits(intr_by_region,selected_traits)

intr_alz@elementMetadata$GENE
intr_alz@elementMetadata$SNPS

intr_alz_by_region@elementMetadata$GENE
intr_alz_by_region@elementMetadata$SNPS

sort(table(getTraits(intr_by_region)), decreasing=TRUE)[1:10]

write.table(sort(table(getTraits(intr_by_region)), decreasing=TRUE),"./results/gwas_traits_by_region.txt")
write.table(intr_by_region,"./results/gwas_intersection_details_by_region.txt")