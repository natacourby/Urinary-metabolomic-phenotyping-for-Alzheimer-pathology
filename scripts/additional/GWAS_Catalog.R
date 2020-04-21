################################################################################
# GWAS CATALOG
################################################################################

#########################################
# Load mQTL results
#########################################
setwd("/home/natalja/Documents/EMIF/KCL")
mets_selected <- read.table("/home/natalja/Documents/EMIF/list_of_important_metabolites.txt",sep="\t",header = T,stringsAsFactors = F)

base <- "UHPOS"
res <- read.table(paste("result_",base,"_annot_full.txt",sep=""),header=T,sep="\t")

res_all <-res[0,]
res_all$base <= ""
for(base in c("UHPOS","URPOS","URNEG")){
  res <- read.table(paste("result_",base,"_annot_full.txt",sep=""),header=T,sep="\t")
  res$SNP <- gsub('_2','',res$SNP)
  res$SNP <- gsub('_1','',res$SNP)
  res$base <- base
  res_all <- rbind(res_all,res)
}

base <-"NMR_Nosey_Urine"
res_nmr <- read.table(paste("result_",base,"_annot_full.txt",sep=""),sep="\t",header=T)
res_nmr$base <- "NMR NOESY"
res_nmr$SNP <- gsub('_2','',res_nmr$SNP)
res_nmr$SNP <- gsub('_1','',res_nmr$SNP)
res_all <- rbind(res_all,res_nmr)

res_genes <- res_all[!is.na(res_all$chr_name) & res_all$chr_name!="",]
res_genes<-res_genes[-grep("CHR_", res_genes$chr_name),]

#########################################
# Load GWAS Catalog from downloaded file 
#########################################
library(gwascat)
gwas <- read.csv("/home/natalja/Documents/R_projects/mQTL/data/gwas_catalog_v1.0-associations_e91_r2018-03-13.tsv",sep="\t")
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

#########################################
# Direct mapping by SNPs
#########################################
intr = ebicat38[ intersect(getRsids(ebicat38), res_genes$SNP) ]
sort(table(getTraits(intr)), decreasing=TRUE)

subsetByTraits(intr,"Schizophrenia")@elementMetadata$SNPS
subsetByTraits(intr,"Amyotrophic lateral sclerosis")@elementMetadata$SNPS


met_list <- 
  unique(res_genes[res_genes$SNP %in% subsetByTraits(intr,"Amyotrophic lateral sclerosis")@elementMetadata$SNPS,]$gene)

intersect(mets_selected$met,met_list)


write.table(sort(table(getTraits(intr)), decreasing=TRUE),"mapping_with_gwas_by_snps.txt",sep="\t")
write.table(intr,"mapping_with_gwas_by_snps2.txt")

#########################################
# Mapping by genomic regions
#########################################
gr6.0 = GRanges(seqnames=res_genes$chr_name,
                IRanges(res_genes$chrom_start, width=1))
mcols(gr6.0)$rsid = as.character(res_genes$SNP)
#seqlevels(gr6.0) = paste("chr", seqlevels(gr6.0), sep="")
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
intr2 = ebicat38[augrs]


intr_alz <- subsetByTraits(intr,selected_traits)

intr2_alz <- subsetByTraits(intr2,selected_traits)

intr_alz@elementMetadata$GENE
intr_alz@elementMetadata$SNPS

intr2_alz@elementMetadata$GENE
intr2_alz@elementMetadata$SNPS

sort(table(getTraits(intr2)), decreasing=TRUE)[1:10]
write.table(sort(table(getTraits(intr2)), decreasing=TRUE),"mapping_with_gwas_by_regions.txt")
write.table(intr2,"mapping_with_gwas_by_regions2.txt")

#########################################
# Plots
#########################################

traitsManh(intr2,
           selr = GRanges(seqnames=res_genes$chr_name,
                          IRanges(res_genes$chrom_start, width=10000000)))

pdf("traits_chr19.pdf")
traitsManh(intr2,
           selr = GRanges(seqnames=4,
                          IRanges(res_genes$chrom_start, width=10000000)),
           traits=getTraits(intr2)
)

dev.off()

traitsManh(intr2,
           selr = GRanges(seqnames=res_genes$chr_name,
                          IRanges(res_genes$chrom_start, width=10000000)),
           traits = getTraits(intr2))


seqlevels(intr2) = paste("chr", seqlevels(intr2), sep="")
seqlevels(alz_traits) = paste("chr", seqlevels(alz_traits), sep="")
seqlevels(gr) = paste("chr", seqlevels(gr), sep="")


# ALMS1 (2: 73,385,758-73,610,793) 
# SLC17A3 (6: 25,833,066-25,882,286) 
# AOX1 (2: 200,585,868-200,677,064) 
# BST1 (4: 15,702,950-15,738,313)
bile_acids <- res_genes[res_genes$gene %in% c("6.13_600.2850m/z","6.13_600.4780m/z","6.13_668.2724m/z","6.13_685.2611m/z","6.13_736.2562m/z","6.30_602.2977m/z","8.01_639.3614n"),]

#parecetamol and diltiazem
paracetaml_diltiazem <- res_genes[res_genes$gene %in% c("2.58_151.0635n","5.01_535.1753m/z", "5.51_551.1705m/z"),]

acisoga <- res_genes[res_genes$gene %in% c("2.62_125.0837n","2.62_92.1133n"),]

butyryl <- res_genes[res_genes$gene %in% c("4.55_190.2275n"),]
methylcytidine_5 <- res_genes[res_genes$gene %in% c("3.80_258.1087m/z","3.80_126.0658m/z"),]
methionine <- res_genes[res_genes$gene %in% c("3.74_249.0860n"),]

methylcytidine_2 <- res_genes[res_genes$gene %in% c("2.03_112.0501m/z"),]

quinine<- res_genes[res_genes$gene %in% c("3.07_283.1580n"),]

amino_isobutyrate  <- res_genes[res_genes$gene %in% c("X1.2051","X1.1952"),]


annotated_met <- res_genes[res_genes$gene %in% c("6.13_600.2850m/z","6.13_600.4780m/z","6.13_668.2724m/z","6.13_685.2611m/z","6.13_736.2562m/z","6.30_602.2977m/z","8.01_639.3614n",
                                                 "2.58_151.0635n","5.01_535.1753m/z", "5.51_551.1705m/z",
                                                 "2.62_125.0837n","2.62_92.1133n",
                                                 "4.55_190.2275n",
                                                 "3.80_258.1087m/z","3.80_126.0658m/z","2.03_112.0501m/z",
                                                 "3.74_249.0860n",
                                                 "2.03_112.0501m/z"),]


  
traitsManh(intr2,
           selr = GRanges(seqnames="chr5",
                          IRanges(amino_isobutyrate$chrom_start, width=10000000)),
           traits=getTraits(ebicat38))


# regions 
# chr19 (0.58e7,0.1e6); 
# chr10 (4.85e7,0.5e6) 
# "Inflammatory skin disease, Crohn's disease, Dietary macronutrient intake,
# Folate pathway vitamin levels,Liver enzyme levels (alkaline phosphatase),
# Multiple sclerosis, Psoriasis, Pubertal anthropometrics, Serum lipase activity,
# Urinary metabolites (H-NMR features), Vitamin B levels in ischemic stroke"                     
# chr19 (1.58e7,0.2e6) Circulating phylloquinone levels
# chr19 (3.27e7, 0.3e6)
# chr4 (6.9e7,0.2e6) obesity-related trait
# chr4 (1.57e7,0.1e6) Parkinson's disease
# chr4 (2.28e7,0.1e6) Blood metabolite levels
# chr14 (5.287e7,0.2e6) Alzheimer's disease (late onset)
# chr14 (7.32e7,0.35e6) Biopolar disoder, Blood metabolites
# chr14 (6.55e7,0.1e6) Menarche
# chr8 (1.84e7,0.1e6) Iron, insulin, metabolic
# chr9 (13.34e7,0.1e6) Liver enzyme levels

chr="chr5"
chr_nr = 5
start_val = 3.498e7
width_val = 0.1e6

onset 
# Colours for Traits
library(RColorBrewer) 
library(ggplot2)
# colors in ggplot2 style
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
cols <- gg_color_hue(length(unique(getTraits(intr2))))

col_vals <- cols[as.factor(getTraits(intr2))]


library(rtracklayer)
#c19tg = import("/home/natalja/Documents/EMIF/Homo_sapiens.GRCh38.90.chromosome.19.gff3")
library(Gviz)
gff3 <- import.gff3("/home/natalja/Documents/EMIF/gencode.v27.basic.annotation.gff3.gz")
gff3_pc <- gff3[gff3$gene_type=="protein_coding",]

gene2symbol <- mcols(gff3_pc)[,c("gene_id","gene_name")]
gene2symbol <- unique(gene2symbol)
rownames(gene2symbol) <- gene2symbol$gene_id
#export.gff3(gff3_pc, "/home/natalja/Documents/EMIF/gencode.v27.basic.annotation_pc.gff3")
txdb <- makeTxDbFromGFF("/home/natalja/Documents/EMIF/gencode.v27.basic.annotation_pc.gff3")




g2 = GRanges(seqnames=chr, IRanges(start=start_val, width=width_val))
#seqlevelsStyle(ebicat38) = "UCSC"


basic = gwcex2gviz_custom3(gene2symbol,txdb,basegr = intr2, contextGR=g2, plot.it=FALSE)



res_genes_19 <- res_genes[res_genes$chr_name==chr_nr,]
res_genes_19 <- amino_isobutyrate[amino_isobutyrate$chr_name==chr_nr,]
mqtl_df<- data.frame(SNP = res_genes_19$SNP,chrom_start = res_genes_19$chrom_start,
                     P_VALUE = -log(res_genes_19$FDR))

mqtl_df2 <- aggregate(mqtl_df,by=list(mqtl_df$SNP),FUN=mean)
mqtl_df2 <- mqtl_df2[,-2]
gr = GRanges(seqnames=chr_nr,
             IRanges(mqtl_df2$chrom_start, width=1))
mcols(gr)$SNPS = as.character(mqtl_df2$Group.1)
values(gr) = mqtl_df2$P_VALUE
genome(gr) <- "hg19"

#gr1 <- gr[gr %in% g2]
#c19ts = c19tg[ which(mcols(c19tg)$type == "Stranger_eqtl") ]
#eqloc = AnnotationTrack(selr, chrom="chr19", genome="hg19", name="Metabolic QTL")

#gr <- GRanges(seqnames="chr7", range=IRanges(start=df$start, end=df$end), strand=str)
mqtlloc <- DataTrack(range=gr, name="Metabolic QTL", genome="hg19")

#gwas_track <- DataTrack(range=range(intr2), name = "GWAS Catalogue", genome="hg19")

displayPars(mqtlloc)$col = "black"
integ = list(basic[[1]],mqtlloc, basic[[2]], basic[[3]])
#width = 4, height = 4, units = 'in', 
png("3-amino-isobutyrate_genomic_level_chr5.png", width = 5, height = 5, units = 'in', res = 300)
plotTracks(integ,from = start_val ,to = start_val+width_val, chromosome = chr)
dev.off()



################################################################################
gwcex2gviz_custom2 = function(gene2symbol,txdb,basegr, 
                              contextGR = GRanges(seqnames="chr17", IRanges(start=37500000, width=1e6)), 
                              txrefpk = "TxDb.Hsapiens.UCSC.hg19.knownGene", genome="hg19",
                              genesympk = "Homo.sapiens",
                              plot.it=TRUE, maxmlp=25 ) {
  #
  # objective is to visualize features of GWAS in gene/transcript context
  #
  # require(Gviz, quietly=TRUE)
  requireNamespace(txrefpk, quietly=TRUE)
  requireNamespace(genesympk, quietly=TRUE)
  # symmap = get(gsub(".db", "SYMBOL", genesympk))
  chrmin = as.character(seqnames(contextGR))
  #
  # the get() here is a hack.  need to have a way of getting relevant object
  # name from txrefpk, indeed allowing more general spec of an exon range resource
  #
  ex = exons( get(txrefpk), columns = c("gene_id", "tx_id", "exon_id"),
              filter=list(exon_chrom = chrmin) )
  txin = ex[ which(overlapsAny(ex, contextGR)) ]
  if (length(txin) == 0) stop("no transcripts in contextGR")
  v = mcols(txin)
  e = v$exon_id
  txl = v$tx_id
  texx = sapply(as.list(v$tx_id), "[", 1)
  g = sapply(as.list(v$gene_id), "[", 1)
  g = unlist(g)
  drop = which(is.na(g))
  k = GRanges(seqnames=chrmin, ranges=ranges(txin), gene=g, exon=e,
              transcript=texx, id=1:length(g))
  if (length(drop) > 0) k = k[-drop]
  kk = mapIds(get(genesympk), keys=mcols(k)$gene, keytype="ENTREZID",
              column="SYMBOL")  # multiVals == first
  mcols(k)$symbol = kk
  GR = GeneRegionTrack(k, chromosome=chrmin, genome=genome)
  
  GR = GeneRegionTrack(txdb, chromosome=chrmin)
  ranges(GR)$symbol <- gene2symbol[ranges(GR)$gene, "gene_name"]
  
  studs = basegr[ which(overlapsAny(basegr, contextGR)) ] 
  traits = unique(getTraits(studs))
  print(traits)
  mlp = mcols(studs)$PVALUE_MLOG
  mlp = ifelse(mlp > maxmlp, maxmlp, mlp)
  mcols(studs)$PVALUE_MLOG = mlp
  sss = mcols(studs)$PVALUE_MLOG
  
  sss2 <- t(matrix(rep(sss,length(traits)),nrow=length(sss)))
  
  i <- 1
  for(trait in traits){
    sss2[i,] <- ifelse(getTraits(studs) != traits[i],NA,sss2[i,])
    i <- i+1
  }
  
  studp = DataTrack( as(studs, "GRanges"), data=sss2, 
                     chromosome=chrmin, genome=genome, name="GWAS Catalog",groups=traits)
  #
  # need to move these controls up to interface
  #
  displayPars(GR)$collapseTranscripts = TRUE
  displayPars(GR)$showId = TRUE
  GR@name = "Genes"
  if (plot.it) plotTracks(list(  studp, GenomeAxisTrack(),  GR))
  invisible( list( studp, GenomeAxisTrack(), GR ) )
}
################################################################################
gwcex2gviz_custom3 = function(gene2symbol,txdb,basegr, 
                              contextGR = GRanges(seqnames="chr17", IRanges(start=37500000, width=1e6)), 
                              txrefpk = "TxDb.Hsapiens.UCSC.hg19.knownGene", genome="hg19",
                              genesympk = "Homo.sapiens",
                              plot.it=TRUE, maxmlp=25, track_name="GWAS Catalog") {
  #
  # objective is to visualize features of GWAS in gene/transcript context
  #
  # require(Gviz, quietly=TRUE)
  requireNamespace(txrefpk, quietly=TRUE)
  requireNamespace(genesympk, quietly=TRUE)
  # symmap = get(gsub(".db", "SYMBOL", genesympk))
  chrmin = as.character(seqnames(contextGR))
  #
  # the get() here is a hack.  need to have a way of getting relevant object
  # name from txrefpk, indeed allowing more general spec of an exon range resource
  #
  ex = exons( get(txrefpk), columns = c("gene_id", "tx_id", "exon_id"),
              filter=list(exon_chrom = chrmin) )
  txin = ex[ which(overlapsAny(ex, contextGR)) ]
  if (length(txin) == 0) stop("no transcripts in contextGR")
  v = mcols(txin)
  e = v$exon_id
  txl = v$tx_id
  texx = sapply(as.list(v$tx_id), "[", 1)
  g = sapply(as.list(v$gene_id), "[", 1)
  g = unlist(g)
  drop = which(is.na(g))
  k = GRanges(seqnames=chrmin, ranges=ranges(txin), gene=g, exon=e,
              transcript=texx, id=1:length(g))
  if (length(drop) > 0) k = k[-drop]
  kk = mapIds(get(genesympk), keys=mcols(k)$gene, keytype="ENTREZID",
              column="SYMBOL")  # multiVals == first
  mcols(k)$symbol = kk
  GR = GeneRegionTrack(k, chromosome=chrmin, genome=genome)
  
  GR = GeneRegionTrack(txdb, chromosome=chrmin)
  ranges(GR)$symbol <- gene2symbol[ranges(GR)$gene, "gene_name"]
  
  studs = basegr[ which(overlapsAny(basegr, contextGR)) ]
  mlp = mcols(studs)$PVALUE_MLOG
  mlp = ifelse(mlp > maxmlp, maxmlp, mlp)
  mcols(studs)$PVALUE_MLOG = mlp
  sss = mcols(studs)$PVALUE_MLOG
  studp = DataTrack( as(studs, "GRanges"), data=sss, 
                     chromosome=chrmin, genome=genome, name=track_name)
  #
  # need to move these controls up to interface
  #
  displayPars(GR)$collapseTranscripts = TRUE
  displayPars(GR)$showId = TRUE
  GR@name = "Genes"
  if (plot.it) plotTracks(list(  studp, GenomeAxisTrack(),  GR))
  invisible( list( studp, GenomeAxisTrack(), GR ) )
}
