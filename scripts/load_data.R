################################################################################
# Script to load selected from metabolic QTL results feature sets and classifiers combinations
# and prepare it for RF algorithm 
# A. Metabolites only;
# B. Metabolites and SNPs;
# C. Metabolites, SNPs and covariates;
# D. Metabolites and covariates.
# Classifier I: 4 original classes: ADC, CTL , sMCI and cMCI; 
# Classifier II: 2 classes when sMCI class is remapped to CTL and cMCI is remapped to ADC;
# Classifier III: 2 classes when both cMCI and sMCI are remapped to ADC;
# Classifier IV: 2 classes when MCIs classes are removed from the analysis.
################################################################################

expr1 <- read.table("/nfs/ma/home/natalja/EMIF/Multimodality_AD_Project/metabolomicQTL/mQTL_data/final_results_urine/selected_matrices/UHPOS_met_sig2.txt",sep="\t",header=T,row.names=1,stringsAsFactors = F)
expr2 <- read.table("/nfs/ma/home/natalja/EMIF/Multimodality_AD_Project/metabolomicQTL/mQTL_data/final_results_urine/selected_matrices/URPOS_met_sig2.txt",sep="\t",header=T,row.names=1,stringsAsFactors = F)
expr3 <- read.table("/nfs/ma/home/natalja/EMIF/Multimodality_AD_Project/metabolomicQTL/mQTL_data/final_results_urine/selected_matrices/URNEG_met_sig2.txt",sep="\t",header=T,row.names=1,stringsAsFactors = F)
expr_nmr <- read.table("/nfs/ma/home/natalja/EMIF/Multimodality_AD_Project/metabolomicQTL/mQTL_data/final_results_urine/selected_matrices/NMR_met_sig2.txt",sep="\t",header=T,row.names=1,stringsAsFactors = F)

expr_nmr <- expr_nmr[,colnames(expr_nmr) %in% colnames(expr1)]
expr1 <- expr1[,colnames(expr1) %in% colnames(expr_nmr)]
expr2 <- expr2[,colnames(expr2) %in% colnames(expr_nmr)]
expr3 <- expr3[,colnames(expr3) %in% colnames(expr_nmr)]

expr_nmr <- expr_nmr[,colnames(expr1)]
expr2 <- expr2[,colnames(expr1)]
expr3 <- expr3[,colnames(expr1)]

expr <- rbind(expr_nmr,expr1)
expr <- rbind(expr,expr2)
expr <- rbind(expr,expr3)

cvrt = read.table("/nfs/ma/home/natalja/EMIF/Multimodality_AD_Project/metabolomicQTL/mQTL_data/covariates/Covariates_all_vertical_text2.txt",sep="\t",header=T,row.names=1)
cvrt <- cvrt[rownames(cvrt) %in% colnames(expr),]
cvrt <- cvrt[colnames(expr),]

################################################################################
####### ORIGINAL DATA ########
expr_orig1 <- read.table(paste("/nfs/ma/home/natalja/EMIF/Multimodality_AD_Project/metabolomicQTL/mQTL_data/results/normalization_final/full_met/UHPOS_original.csv",sep=""),sep=",",header=T,row.names=1,stringsAsFactors = F)
epxr_orig1 <- expr_orig1[rownames(expr_orig1) %in% rownames(expr1),]
expr_orig2 <- read.table(paste("/nfs/ma/home/natalja/EMIF/Multimodality_AD_Project/metabolomicQTL/mQTL_data/results/normalization_final/full_met/URPOS_original.csv",sep=""),sep=",",header=T,row.names=1,stringsAsFactors = F)
epxr_orig2 <- expr_orig2[rownames(expr_orig2) %in% rownames(expr2),]
expr_orig3 <- read.table(paste("/nfs/ma/home/natalja/EMIF/Multimodality_AD_Project/metabolomicQTL/mQTL_data/results/normalization_final/full_met/URNEG_original.csv",sep=""),sep=",",header=T,row.names=1,stringsAsFactors = F)
epxr_orig3 <- expr_orig3[rownames(expr_orig3) %in% rownames(expr3),]
expr_orig_nmr <- read.table("/nfs/ma/home/natalja/EMIF/Multimodality_AD_Project/metabolomicQTL/mQTL_data/results/normalization_final/full_met/NMR_Nosey_Urine_original.txt",sep="\t",header=T,row.names=1,stringsAsFactors = F)

expr_orig1 <- expr_orig1[rownames(expr_orig1) %in% rownames(expr1),]
expr_orig2 <- expr_orig2[rownames(expr_orig2) %in% rownames(expr2),]
expr_orig3 <- expr_orig3[rownames(expr_orig3) %in% rownames(expr3),]
expr_orig_nmr <- expr_orig_nmr[rownames(expr_orig_nmr) %in% rownames(expr_nmr),]

expr_orig_nmr <- expr_orig_nmr[,colnames(expr_orig_nmr) %in% colnames(expr_orig1)]
expr_orig1 <- expr_orig1[,colnames(expr_orig1) %in% colnames(expr_orig_nmr)]
expr_orig2 <- expr_orig2[,colnames(expr_orig2) %in% colnames(expr_orig_nmr)]
expr_orig3 <- expr_orig3[,colnames(expr_orig3) %in% colnames(expr_orig_nmr)]

expr_orig_nmr <- expr_orig_nmr[,colnames(expr_orig1)]
expr_orig2 <- expr_orig2[,colnames(expr_orig1)]
expr_orig3 <- expr_orig3[,colnames(expr_orig1)]

expr_orig <- rbind(expr_orig_nmr,expr_orig1)
expr_orig <- rbind(expr_orig,expr_orig2)
expr_orig <- rbind(expr_orig,expr_orig3)

expr_orig <- expr_orig[,colnames(expr_orig) %in% colnames(expr)]
expr_orig <- expr_orig[,colnames(expr)]

QN_normalize <- function(met){
  for (i in 1:nrow(met)){
    mat <- met[i,]
    mat = t(apply(mat, 1, rank, ties.method = "average"));
    mat = qnorm(mat / (ncol(met)+1))
    met[i,] <- mat 
  }
  return(met)
}

expr_orig_qn <- QN_normalize(expr_orig)
#####################################################################
#####################################################################


# a) Metabolites only and d) Metabolites and covariates
df <- as.data.frame(t(expr))
df$Gender <- cvrt$Gender
df$Centre <- cvrt$Centre
df$Age <- cvrt$Age
df$Diagnosis1 <- as.character(cvrt$Diagnosis)
df$Diagnosis2 <- as.character(cvrt$Diagnosis)
df$Diagnosis3 <- as.character(cvrt$Diagnosis)
df$Diagnosis2[df$Diagnosis2=="cMCI"] <- "ADC"
df$Diagnosis2[df$Diagnosis2=="sMCI"] <- "CTL"
df$Diagnosis3[df$Diagnosis3=="cMCI"] <- "ADC"
df$Diagnosis3[df$Diagnosis3=="sMCI"] <- "ADC"

# a) 
drops <- c("Gender","Centre","Age","Diagnosis2","Diagnosis3")
df1_a <- df[,!(names(df) %in% drops)]
names(df1_a)[names(df1_a) == 'Diagnosis1'] <- 'Diagnosis'

drops <- c("Gender","Centre","Age","Diagnosis1","Diagnosis3")
df2_a <- df[,!(names(df) %in% drops)]
names(df2_a)[names(df2_a) == 'Diagnosis2'] <- 'Diagnosis'

drops <- c("Gender","Centre","Age","Diagnosis1","Diagnosis2")
df3_a <- df[,!(names(df) %in% drops)]
names(df3_a)[names(df3_a) == 'Diagnosis3'] <- 'Diagnosis'

drops <- c("Gender","Centre","Age","Diagnosis2","Diagnosis3")
df4_a <- df[cvrt$Diagnosis %in% c("CTL","ADC"),!(names(df) %in% drops)]
names(df4_a)[names(df4_a) == 'Diagnosis1'] <- 'Diagnosis'

# or set d
drops <- c("Diagnosis2","Diagnosis3")
df1_d <- df[,!(names(df) %in% drops)]
names(df1_d)[names(df1_d) == 'Diagnosis1'] <- 'Diagnosis'

drops <- c("Diagnosis1","Diagnosis3")
df2_d <- df[,!(names(df) %in% drops)]
names(df2_d)[names(df2_d) == 'Diagnosis2'] <- 'Diagnosis'

drops <- c("Diagnosis1","Diagnosis2")
df3_d <- df[,!(names(df) %in% drops)]
names(df3_d)[names(df3_d) == 'Diagnosis3'] <- 'Diagnosis'

drops <- c("Diagnosis2","Diagnosis3")
df4_d <- df[cvrt$Diagnosis %in% c("CTL","ADC"),!(names(df) %in% drops)]
names(df4_d)[names(df4_d) == 'Diagnosis1'] <- 'Diagnosis'

# set b
dna_all <- read.table("/nfs/ma/home/natalja/EMIF/Multimodality_AD_Project/metabolomicQTL/mQTL_data/final_results_urine/selected_matrices/dna_all.txt",sep="\t",header=T,row.names=1,stringsAsFactors = F)
# The common samples for NMR and LC-MS
dna_all<- dna_all[,colnames(dna_all) %in% colnames(expr)]
expr_dna <- expr[,colnames(expr) %in% colnames(dna_all)]
dna_all <- dna_all[,colnames(expr_dna)]
dna_all [is.na(dna_all )] <- 0
cvrt = read.table("/nfs/ma/home/natalja/EMIF/Multimodality_AD_Project/metabolomicQTL/mQTL_data/covariates/Covariates_all_vertical_text2.txt",sep="\t",header=T,row.names=1)
cvrt <- cvrt[rownames(cvrt) %in% colnames(expr_dna),]
cvrt <- cvrt[colnames(expr_dna),]

df <- as.data.frame(cbind(t(expr_dna),t(dna_all)))
df$Gender <- cvrt$Gender
df$Centre <- cvrt$Centre
df$Age <- cvrt$Age
df$Diagnosis1 <- as.character(cvrt$Diagnosis)
df$Diagnosis2 <- as.character(cvrt$Diagnosis)
df$Diagnosis3 <- as.character(cvrt$Diagnosis)
df$Diagnosis2[df$Diagnosis2=="cMCI"] <- "ADC"
df$Diagnosis2[df$Diagnosis2=="sMCI"] <- "CTL"
df$Diagnosis3[df$Diagnosis3=="cMCI"] <- "ADC"
df$Diagnosis3[df$Diagnosis3=="sMCI"] <- "ADC"

# b) 
drops <- c("Gender","Centre","Age","Diagnosis2","Diagnosis3")
df1_b <- df[,!(names(df) %in% drops)]
names(df1_b)[names(df1_b) == 'Diagnosis1'] <- 'Diagnosis'

drops <- c("Gender","Centre","Age","Diagnosis1","Diagnosis3")
df2_b <- df[,!(names(df) %in% drops)]
names(df2_b)[names(df2_b) == 'Diagnosis2'] <- 'Diagnosis'

drops <- c("Gender","Centre","Age","Diagnosis1","Diagnosis2")
df3_b <- df[,!(names(df) %in% drops)]
names(df3_b)[names(df3_b) == 'Diagnosis3'] <- 'Diagnosis'

drops <- c("Gender","Centre","Age","Diagnosis2","Diagnosis3")
df4_b <- df[cvrt$Diagnosis %in% c("CTL","ADC"),!(names(df) %in% drops)]
names(df4_b)[names(df4_b) == 'Diagnosis1'] <- 'Diagnosis'

# or set c
drops <- c("Diagnosis2","Diagnosis3")
df1_c <- df[,!(names(df) %in% drops)]
names(df1_c)[names(df1_c) == 'Diagnosis1'] <- 'Diagnosis'

drops <- c("Diagnosis1","Diagnosis3")
df2_c <- df[,!(names(df) %in% drops)]
names(df2_c)[names(df2_c) == 'Diagnosis2'] <- 'Diagnosis'

drops <- c("Diagnosis1","Diagnosis2")
df3_c <- df[,!(names(df) %in% drops)]
names(df3_c)[names(df3_c) == 'Diagnosis3'] <- 'Diagnosis'

drops <- c("Diagnosis2","Diagnosis3")
df4_c <- df[cvrt$Diagnosis %in% c("CTL","ADC"),!(names(df) %in% drops)]
names(df4_c)[names(df4_c) == 'Diagnosis1'] <- 'Diagnosis'

df1_a$Diagnosis<- as.factor(df1_a$Diagnosis)
df1_b$Diagnosis<- as.factor(df1_b$Diagnosis)
df1_c$Diagnosis<- as.factor(df1_c$Diagnosis)
df1_d$Diagnosis<- as.factor(df1_d$Diagnosis)
df2_a$Diagnosis<- as.factor(df2_a$Diagnosis)
df2_b$Diagnosis<- as.factor(df2_b$Diagnosis)
df2_c$Diagnosis<- as.factor(df2_c$Diagnosis)
df2_d$Diagnosis<- as.factor(df2_d$Diagnosis)
df3_a$Diagnosis<- as.factor(df3_a$Diagnosis)
df3_b$Diagnosis<- as.factor(df3_b$Diagnosis)
df3_c$Diagnosis<- as.factor(df3_c$Diagnosis)
df3_d$Diagnosis<- as.factor(df3_d$Diagnosis)
df4_a$Diagnosis<- as.factor(df4_a$Diagnosis)
df4_b$Diagnosis<- as.factor(df4_b$Diagnosis)
df4_c$Diagnosis<- as.factor(df4_c$Diagnosis)
df4_d$Diagnosis<- as.factor(df4_d$Diagnosis)
names(df1_a) <- make.names(names(df1_a))
names(df1_b) <- make.names(names(df1_b))
names(df1_c) <- make.names(names(df1_c))
names(df1_d) <- make.names(names(df1_d))
names(df2_a) <- make.names(names(df2_a))
names(df2_b) <- make.names(names(df2_b))
names(df2_c) <- make.names(names(df2_c))
names(df2_d) <- make.names(names(df2_d))
names(df3_a) <- make.names(names(df3_a))
names(df3_b) <- make.names(names(df3_b))
names(df3_c) <- make.names(names(df3_c))
names(df3_d) <- make.names(names(df3_d))
names(df4_a) <- make.names(names(df4_a))
names(df4_b) <- make.names(names(df4_b))
names(df4_c) <- make.names(names(df4_c))
names(df4_d) <- make.names(names(df4_d)) 


