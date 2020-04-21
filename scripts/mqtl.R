mqtl_analysis <- function(snps,met_matrix_name,cov_name,output_file_name){ 
  require(MatrixEQTL)
  #base.dir = "/home/natalja/Documents/EMIF/KCL/mQTL/"
  #setwd(base.dir)
  
  useModel = modelLINEAR
  #SNP_file_name = paste(base.dir, dna_matrix_name, sep="")
  #snps_location_file_name = paste(base.dir, "/dna_matrix_pos.tsv", sep="");
  metabolomics_file_name = met_matrix_name #paste(base.dir, met_matrix_name, sep="")
  #gene_location_file_name = paste(base.dir, "/metloc.txt", sep="");
  covariates_file_name = cov_name #paste(base.dir, cov_name, sep="")
  
  pvOutputThreshold = 1e-5;
  errorCovariance = numeric();
  min.pv.by.genesnp = TRUE 
  
  tic_load = proc.time()[3];
  
  #snps = SlicedData$new();
  #snps$fileDelimiter = "\t"; # the TAB character
  #snps$fileOmitCharacters = 'NA' ;# denote missing values;
  #snps$fileSkipRows = 1; # one row of column labels
  #snps$fileSkipColumns = 1; # one column of row labels
  #snps$fileSliceSize = 10000; # read file in pieces of 10,000 rows
  #snps$LoadFile(SNP_file_name);
  
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
  #cat('eQTL time: ', toc_load-tic_load, ' sec\n');
  
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

# Auxiliary functions
##################################################################################################
extract_features <- function(file_name, new_name){
  met_data <- read.table(file_name,header = TRUE,stringsAsFactors = FALSE)
  features <- read.table("/home/natalja/Documents/EMIF/KCL/mQTL/mQTL_data/features.txt",stringsAsFactors = FALSE)
  features <- as.vector(features$V1)
  features <- as.character(features)
  if(length(intersect(rownames(met_data),features))==0){
    rownames(met_data)<-met_data[,1]
  }
  else {
    features <- intersect(rownames(met_data),features)
  }
  met_data <- met_data[features,]
  write.table(met_data,new_name,sep="\t",quote = FALSE)
}
########################
extract_features_results <- function(file_name, new_name){
  met_data <- read.table(file_name,header = TRUE,stringsAsFactors = FALSE)
  features <- read.table("/home/natalja/Documents/EMIF/KCL/mQTL/mQTL_data/features.txt",stringsAsFactors = FALSE)
  features <- as.vector(features$V1)
  features <- as.character(features)
  features <- intersect(met_data$metabolite,features)
  if (length(features)>0){
    met_data <- met_data[met_data$metabolite %in% features,]
    write.table(met_data,new_name,sep="\t",quote = FALSE)
  }
  else{
    print("No such features!")
  }
}
########################
match_matrices <- function(matrix1_name, matrix2_name, matrix2_new_name){
  matrix1 <- read.table(matrix1_name, header = TRUE,stringsAsFactors = FALSE ,skip=0, nrows=1)
  matrix2 <- read.table(matrix2_name,header = TRUE,stringsAsFactors = FALSE)
  matrix2 <- matrix2[,intersect(colnames(matrix1),colnames(matrix2))]
  samples <- colnames(matrix2) 
  matrix2 <- cbind(rownames(matrix2),matrix2)
  colnames(matrix2) <- c("id",samples)
  write.table(matrix2,matrix2_new_name,sep="\t",quote = FALSE,row.names = FALSE)
}
########################
check_matrices <- function(dna_matrix_name,met_matrix_name,cov_name){
  matrix1 <- read.table(dna_matrix_name, header = TRUE,stringsAsFactors = FALSE ,skip=0, nrows=1)
  matrix2 <- read.table(met_matrix_name, header = TRUE,stringsAsFactors = FALSE ,skip=0, nrows=1)
  matrix3 <- read.table(cov_name, header = TRUE,stringsAsFactors = FALSE ,skip=0, nrows=1)
  a <- all.equal(colnames(matrix1)[-1],colnames(matrix2)[-1])
  b <- all.equal(colnames(matrix2)[-1],colnames(matrix3)[-1])
  c <- all.equal(colnames(matrix1)[-1],colnames(matrix3)[-1])
  return(c(a,b,c))
}
########################
# Bonferroni correction
bonferroni_correction <- function(file_name){
  mqtl <- read.table(file_name,header = TRUE,stringsAsFactors = FALSE)
  dim(mqtl)
  mqtl$bonf <- p.adjust(mqtl$p.value,"bonferroni")
  mqtl <- mqtl[mqtl$bonf<0.01,]
  dim(mqtl)
  names(mqtl) <- c("SNP","metabolite","beta", "t.stat","p.value" ,"FDR","bonf") 
  mqtl <- mqtl[-c(6)]
  write.table(mqtl,file_name,sep="\t",row.names = FALSE)
}
########################
# Permutation
permute_matrix <- function(file_input_name,file_output_name){
  met_data <- read.table(file_input_name,header = TRUE,stringsAsFactors = FALSE)
  met_data2 <- met_data[,c(1,sample(2:ncol(met_data)))]
  colnames(met_data2) <- colnames(met_data)
  write.table(met_data2,file_output_name,sep="\t",row.names = FALSE)
}
##################################################################################################

# Manhattam plots
library(lattice)
manhattan.plot<-function(chr, pos, pvalue, 
                         sig.level=NA, annotate=NULL, ann.default=list(),
                         should.thin=T, thin.pos.places=2, thin.logp.places=2, 
                         xlab="Chromosome", ylab=expression(-log[10](p-value)),
                         col=c("gray","darkgray"), panel.extra=NULL, pch=20, cex=0.8,...) {
  
  if (length(chr)==0) stop("chromosome vector is empty")
  if (length(pos)==0) stop("position vector is empty")
  if (length(pvalue)==0) stop("pvalue vector is empty")
  
  #make sure we have an ordered factor
  if(!is.ordered(chr)) {
    chr <- ordered(chr)
  } else {
    chr <- chr[,drop=T]
  }
  
  #make sure positions are in kbp
  if (any(pos>1e6)) pos<-pos/1e6;
  
  #calculate absolute genomic position
  #from relative chromosomal positions
  posmin <- tapply(pos,chr, min);
  posmax <- tapply(pos,chr, max);
  posshift <- head(c(0,cumsum(posmax)),-1);
  names(posshift) <- levels(chr)
  genpos <- pos + posshift[chr];
  getGenPos<-function(cchr, cpos) {
    p<-posshift[as.character(cchr)]+cpos
    return(p)
  }
  
  #parse annotations
  grp <- NULL
  ann.settings <- list()
  label.default<-list(x="peak",y="peak",adj=NULL, pos=3, offset=0.5, 
                      col=NULL, fontface=NULL, fontsize=NULL, show=F)
  parse.label<-function(rawval, groupname) {
    r<-list(text=groupname)
    if(is.logical(rawval)) {
      if(!rawval) {r$show <- F}
    } else if (is.character(rawval) || is.expression(rawval)) {
      if(nchar(rawval)>=1) {
        r$text <- rawval
      }
    } else if (is.list(rawval)) {
      r <- modifyList(r, rawval)
    }
    return(r)
  }
  
  if(!is.null(annotate)) {
    if (is.list(annotate)) {
      grp <- annotate[[1]]
    } else {
      grp <- annotate
    } 
    if (!is.factor(grp)) {
      grp <- factor(grp)
    }
  } else {
    grp <- factor(rep(1, times=length(pvalue)))
  }
  
  ann.settings<-vector("list", length(levels(grp)))
  ann.settings[[1]]<-list(pch=pch, col=col, cex=cex, fill=col, label=label.default)
  
  if (length(ann.settings)>1) { 
    lcols<-trellis.par.get("superpose.symbol")$col 
    lfills<-trellis.par.get("superpose.symbol")$fill
    for(i in 2:length(levels(grp))) {
      ann.settings[[i]]<-list(pch=pch, 
                              col=lcols[(i-2) %% length(lcols) +1 ], 
                              fill=lfills[(i-2) %% length(lfills) +1 ], 
                              cex=cex, label=label.default);
      ann.settings[[i]]$label$show <- T
    }
    names(ann.settings)<-levels(grp)
  }
  for(i in 1:length(ann.settings)) {
    if (i>1) {ann.settings[[i]] <- modifyList(ann.settings[[i]], ann.default)}
    ann.settings[[i]]$label <- modifyList(ann.settings[[i]]$label, 
                                          parse.label(ann.settings[[i]]$label, levels(grp)[i]))
  }
  if(is.list(annotate) && length(annotate)>1) {
    user.cols <- 2:length(annotate)
    ann.cols <- c()
    if(!is.null(names(annotate[-1])) && all(names(annotate[-1])!="")) {
      ann.cols<-match(names(annotate)[-1], names(ann.settings))
    } else {
      ann.cols<-user.cols-1
    }
    for(i in seq_along(user.cols)) {
      if(!is.null(annotate[[user.cols[i]]]$label)) {
        annotate[[user.cols[i]]]$label<-parse.label(annotate[[user.cols[i]]]$label, 
                                                    levels(grp)[ann.cols[i]])
      }
      ann.settings[[ann.cols[i]]]<-modifyList(ann.settings[[ann.cols[i]]], 
                                              annotate[[user.cols[i]]])
    }
  }
  rm(annotate)
  
  #reduce number of points plotted
  if(should.thin) {
    logp <- round(-log10(pvalue),thin.logp.places)
    thinned <- unique(data.frame(
      logp=logp <- logp - min(logp), 
      pos=round(genpos,thin.pos.places), 
      chr=chr,
      grp=grp)
    )
    logp <- thinned$logp
    genpos <- thinned$pos
    chr <- thinned$chr
    grp <- thinned$grp
    rm(thinned)
  } else {
    logp <- -log10(pvalue)
    logp <- logp - min(logp)
  }
  rm(pos, pvalue)
  gc()
  
  #custom axis to print chromosome names
  axis.chr <- function(side,...) {
    if(side=="bottom") {
      panel.axis(side=side, outside=T,
                 at=((posmax+posmin)/2+posshift),
                 labels=levels(chr), 
                 ticks=F, rot=0,
                 check.overlap=F
      )
    } else if (side=="top" || side=="right") {
      panel.axis(side=side, draw.labels=F, ticks=F);
    }
    else {
      axis.default(side=side,...);
    }
  }
  
  #make sure the y-lim covers the range (plus a bit more to look nice)
  prepanel.chr<-function(x,y,...) { 
    A<-list();
    maxy<-ceiling(max(y, ifelse(!is.na(sig.level), -log10(sig.level), 0)))+.5;
    A$ylim=c(0,maxy);
    A;
  }
  
  xyplot(logp~genpos, chr=chr, groups=grp,
         axis=axis.chr, ann.settings=ann.settings, 
         prepanel=prepanel.chr, scales=list(axs="i"),
         panel=function(x, y, ..., getgenpos) {
           if(!is.na(sig.level)) {
             #add significance line (if requested)
             panel.abline(h=-log10(sig.level), lty=2);
           }
           panel.superpose(x, y, ..., getgenpos=getgenpos);
           if(!is.null(panel.extra)) {
             panel.extra(x,y, getgenpos, ...)
           }
         },
         panel.groups = function(x,y,..., subscripts, group.number) {
           A<-list(...)
           #allow for different annotation settings
           gs <- ann.settings[[group.number]]
           A$col.symbol <- gs$col[(as.numeric(chr[subscripts])-1) %% length(gs$col) + 1]    
           A$cex <- gs$cex[(as.numeric(chr[subscripts])-1) %% length(gs$cex) + 1]
           A$pch <- gs$pch[(as.numeric(chr[subscripts])-1) %% length(gs$pch) + 1]
           A$fill <- gs$fill[(as.numeric(chr[subscripts])-1) %% length(gs$fill) + 1]
           A$x <- x
           A$y <- y
           do.call("panel.xyplot", A)
           #draw labels (if requested)
           if(gs$label$show) {
             gt<-gs$label
             names(gt)[which(names(gt)=="text")]<-"labels"
             gt$show<-NULL
             if(is.character(gt$x) | is.character(gt$y)) {
               peak = which.max(y)
               center = mean(range(x))
               if (is.character(gt$x)) {
                 if(gt$x=="peak") {gt$x<-x[peak]}
                 if(gt$x=="center") {gt$x<-center}
               }
               if (is.character(gt$y)) {
                 if(gt$y=="peak") {gt$y<-y[peak]}
               }
             }
             if(is.list(gt$x)) {
               gt$x<-A$getgenpos(gt$x[[1]],gt$x[[2]])
             }
             do.call("panel.text", gt)
           }
         },
         xlab=xlab, ylab=ylab, 
         panel.extra=panel.extra, getgenpos=getGenPos, ...
  );
}

