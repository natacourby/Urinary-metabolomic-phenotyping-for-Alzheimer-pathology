library(biomaRt)
library(MatrixEQTL)

args <- commandArgs(trailingOnly=TRUE)
base <- args[1]
diagnosis <- ""


source("mqtl.R")
# load snps 
SNP_file_name <- paste("./data/dna_matrices/dna_matrix_",base,diagnosis,".tsv",sep="") 
snps = SlicedData$new();
snps$fileDelimiter = "\t"; # the TAB character
snps$fileOmitCharacters = 'NA' ;# denote missing values;
snps$fileSkipRows = 1; # one row of column labels
snps$fileSkipColumns = 1; # one column of row labels
snps$fileSliceSize = 10000; # read file in pieces of 10,000 rows
snps$LoadFile(SNP_file_name);

met_matrix_name <-paste("./data/metabolomics_matrices/normalised/",base,diagnosis,".tsv",sep="")
cov_name <- paste("./data/covariates/Covariates_",base,diagnosis,".txt",sep="") 
output_file_name <- paste("./results/mqtl_",base,diagnosis,".txt",sep="") 
  
# mQTL Analysis
mqtl_analysis(snps,met_matrix_name,cov_name,output_file_name)

# To remove false positives (acquired by chance QTL results), we repeated the analysis on permuted metabolomic data.
permuted_met_matrix_name <- paste(gsub(".tsv","",met_matrix_name),"_permuted.tsv",sep="")
permute_matrix(met_matrix_name,permuted_met_matrix_name)
permuted_output_file_name <- paste("./results/mqtl_permuted_",base,diagnosis,".txt",sep="") 

# mQTL analysis on permuted data
mqtl_analysis(snps,permuted_met_matrix_name,cov_name,permuted_output_file_name)
  
res <- read.table(output_file_name,header=T,sep="\t")
res_permuted <- read.table(permuted_output_file_name,header=T,sep="\t")
res_sig <- res[res$FDR<0.01,]

if (nrow(res_permuted[res_permuted$FDR<0.01,])==0){ # Proceed only if permuted results are not significant# Proceed only if permuted results are not significant
    res_sig$SNP <- gsub('_2','',res_sig$SNP)
    res_sig$SNP <- gsub('_1','',res_sig$SNP)
    
    # Get annotations
    snp.db <- useMart("ENSEMBL_MART_SNP", dataset="hsapiens_snp")
    nt.biomart <- getBM(c("refsnp_id","allele","chr_name","chrom_start","ensembl_gene_stable_id"),
                        filters="snp_filter",
                        values=res_sig$SNP,
                        mart=snp.db)
    genes <- unique(nt.biomart$ensembl_gene_stable_id)
    ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
    genes_annot <- getBM(attributes=c('ensembl_gene_id','hgnc_symbol','chromosome_name','description'), 
                         filters = 'ensembl_gene_id', values = genes, mart = ensembl)
    res_annot <- merge(res_sig,nt.biomart,by.x="SNP",by.y="refsnp_id",all.x=T,all.y=F)
    res_annot <- merge(res_annot,genes_annot,by.x="ensembl_gene_stable_id",by.y="ensembl_gene_id",all.x=T,all.y=F)
    res_annot$nr <- 0
    res <- res_annot[0,]
    for(metabolite in as.character(unique(res_annot$gene))){
      met_group <- res_annot[res_annot$gene==metabolite,]
      for(chr in as.character(unique(met_group$chr_name))){
        if(!(is.na(chr))){
          met_chr_group <-  met_group[met_group$chr_name==chr,]
          snp <- met_chr_group[which.min(met_chr_group$FDR),]
          snp$nr <- nrow(met_chr_group)
          res <- rbind(res,snp)
        }
      }
    }
    # unique per gene/metabolite pair results
    write.table(res,paste("./results/result_",base,diagnosis,"_annot.txt",sep="") ,sep="\t",row.names = FALSE)
    # unique per SNP/metabolite pair results
    write.table(res_annot,paste("./results/result_",base,diagnosis,"_annot_full.txt",sep=""),sep="\t",row.names = FALSE)
}


