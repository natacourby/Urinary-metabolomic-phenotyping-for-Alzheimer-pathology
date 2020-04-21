################################################################################
# Heatmap from selected metabolites
################################################################################
# feature selection by correlation (WEKA result)
library(gplots)
library(RColorBrewer) 
library(ggplot2)
# colors in ggplot2 style
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
hmcol <- brewer.pal(11,"RdBu")
as.fumeric <- function(x,levels=unique(x)) {
  if(!is.character(x)) stop("'x' must be a character")
  as.numeric(factor(x,levels=levels))
}

source("/hps/nobackup/ma/natalja/data/mQTL_data/scripts/load_data.R")



method <- "4D"#"corr" # "LR"
selected_mets <- read.table(paste("/hps/nobackup/ma/natalja/data/mQTL_data/final_results_urine/feature_selection_",method,".txt",sep=""),sep="\t", header=FALSE, stringsAsFactors=F)

selected_mets <- read.table("/hps/nobackup/ma/natalja/data/mQTL_data/final_results_urine/feature_selection_annotated_ANOVA.txt",sep="\t", header=TRUE, stringsAsFactors=F)
expr_selected <- expr[rownames(expr) %in% selected_mets$Feature,]

# Subset from normalized metabolite matrix 
expr_selected <- expr[rownames(expr) %in% selected_mets$V1,]
#expr_selected_nmr <- expr_nmr[rownames(expr_nmr) %in% selected_mets$V1,]
#expr_selected <- rbind(expr_selected_lcms,expr_selected_nmr)
cvrt = read.table("/hps/nobackup/ma/natalja/data/mQTL_data/covariates/Covariates_all_vertical_text2.txt",sep="\t",header=T,row.names=1)
cvrt <- cvrt[rownames(cvrt) %in% colnames(expr_selected),]
cvrt <- cvrt[colnames(expr_selected),]
df <- as.data.frame(t(expr_selected))
df$Diagnosis <- as.character(cvrt$Diagnosis)

cvrt$Diagnosis<-as.character(cvrt$Diagnosis)
cvrt[cvrt$Diagnosis=="ADC",]$Diagnosis <- "AD"
cvrt$Diagnosis<-as.factor(cvrt$Diagnosis)
###############################
# Heatmap of selected features 
###############################

diagnoses <- c("AD", "cMCI", "CTL","sMCI")
cols <- c("#00AFBB", "#E7B800", "#FC4E07","#4F7942")
#re-order levels to correspond to ggplot2 colors and order
cvrt$Diagnosis <- factor(cvrt$Diagnosis, levels = diagnoses)

m <- as.matrix(expr_selected[order(cvrt$Diagnosis)])
rownames(selected_mets) <- selected_mets$Feature
selected_mets <- selected_mets[rownames(selected_mets) %in% rownames(m),]
selected_mets <- selected_mets[rownames(m),]
rownames(m) <- selected_mets$Annotation

png("heatmap_annotated.png", width = 10, height = 9, units = 'in', res=300)
t <- heatmap.2(m,colsep=c(162, 198, 331), sepwidth=2 , sepcolor=c("black","black","black"), col=brewer.pal(11,"RdBu"), 
               trace="none",Colv="none", dendrogram = "row",cexRow=0.9,#key=F, 
               colCol=cols[cvrt[order(cvrt$Diagnosis),]$Diagnosis],margins=c(10,15),key=TRUE, keysize = 1, symkey=FALSE, density.info="none",srtRow=45)

# par(lend = 1)           # square line ends for the color legend
# legend("topleft",      # location of the legend on the heatmap plot
#        legend = diagnoses, # category labels
#        col = cols,  # color key
#        lty= 1,             # line style
#        lwd = 10            # line width
# )
dev.off()
