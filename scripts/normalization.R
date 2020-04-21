# Normalization for MS data

##################################################################################################
# Quantile Normalization of metabolomic data
##################################################################################################
QN_normalize <- function(file_name, output_file_name){
  met <-read.table(file_name,header=T,sep="\t")
  for (i in 1:nrow(met)){
    mat <- met[i,]
    mat = t(apply(mat, 1, rank, ties.method = "average"));
    mat = qnorm(mat / (ncol(met)+1))
    met[i,] <- mat 
  }
  write.table(met,output_file_name,sep="\t")
}

##################################################################################################
# EigenMS Normalization of metabolomic data
##################################################################################################
EigenMS_normalize <- function(file_name, output_file_name, covariates_to_preserve_file_name){
  # source in the EigenMS and correlation heatmap functions
  source("EigenMS.R")
  
  cov <- read.table(covariates_to_preserve_file_name,header = T,sep="\t")
  rownames(cov) <- cov$Sample
  
  ddata = read.table(file_name, header=TRUE, row.names=1, sep=",") 
  #ddata <- ddata[,-(2:5)]
  ddata <- ddata[,intersect(colnames(ddata),rownames(cov))] 
  ddata <- ddata[,order(intersect(colnames(ddata),rownames(cov)))]
  
  cov <- cov[intersect(colnames(ddata),rownames(cov)),]
  cov <- cov[order(intersect(colnames(ddata),rownames(cov))),]
  grps2 <- as.factor(cov$EigenMS_factor)
  
  dim(ddata) # 2802 metabolites x 227 samples
  scaleShift=abs(min(ddata, na.rm = TRUE))+1  
  m_logInts = log(ddata+scaleShift)
  
  m_nummiss = sum(is.na(m_logInts)) # 2000 total values, 700 missing values
  
  m_numtot = dim(m_logInts)[1] * dim(m_logInts)[2] # 8000 total observations
  m_percmiss = m_nummiss/m_numtot  # 8.75% percent missing observations
  
  m_prot.info = cbind(rownames(ddata),rownames(ddata)) # all unique metabolite/peptide IDs, otherwise not possible to have as row names

  # Example part 1 - single factor
  m_ints_eig1 = eig_norm1(m=m_logInts,treatment=grps2,prot.info=m_prot.info,write_to_file = output_file_name)
  m_ints_norm1 = eig_norm2(rv=m_ints_eig1) 
  
  write.table(m_ints_norm1$norm_m,output_file_name,sep="\t")
}
##################################################################################################
##################################################################################################