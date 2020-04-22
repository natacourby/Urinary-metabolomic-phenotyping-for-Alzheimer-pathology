##################################################################################################
# mQTL analysis function
##################################################################################################
mqtl_analysis <- function(snps,met_matrix_name,cov_name,output_file_name){ 
  require(MatrixEQTL)
  useModel = modelLINEAR
  metabolomics_file_name = met_matrix_name
  covariates_file_name = cov_name
  pvOutputThreshold = 1e-5;
  errorCovariance = numeric();
  min.pv.by.genesnp = TRUE 
  tic_load = proc.time()[3];
  
  ## Load metabolomics data
  metabolites = SlicedData$new();
  metabolites$fileDelimiter = '\t'; # the TAB character
  metabolites$fileOmitCharacters = 'NA'; # denote missing values;
  metabolites$fileSkipRows = 1; # one row of column labels
  metabolites$fileSkipColumns = 1; # one column of row labels
  metabolites$fileSliceSize = 10000; # read file in pieces of 10,000 rows
  metabolites$LoadFile(metabolomics_file_name);
  
  ## Load covariates
  cvrt = SlicedData$new();
  cvrt$fileDelimiter = '\t'; # the TAB character
  cvrt$fileOmitCharacters = 'NA'; # denote missing values;
  cvrt$fileSkipRows = 1; # one row of column labels
  cvrt$fileSkipColumns = 1; # one column of row labels
  cvrt$fileSliceSize = snps$nCols()+1; # read file in one piece
  if(length(covariates_file_name)>0) {
    cvrt$LoadFile(covariates_file_name);
  }
  
  toc_load = proc.time()[3];
  
  ## Run the analysis
  {
    tic_mqtl = proc.time()[3];
    me=Matrix_eQTL_engine(snps, metabolites,     cvrt,output_file_name,pvOutputThreshold,useModel, errorCovariance, verbose=TRUE,pvalue.hist = "qqplot");
    # pvalue.hist = 100
    toc_mqtl = proc.time()[3];
  }
  #plot(me, col="grey")
  #plot(me, pch = 16, cex = 0.7)
}  
