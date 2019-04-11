# Note: This Script is mostly based on the instruction manual of RaceID2 and StemID.
install.packages("dplyr")
library(dplyr)

setwd("D:/Userdata/jj.koning/MIKE/lymphnode-sc-transcriptomics-seperatescripts/data")
annotations <- read.csv("annotations.tsv", sep = "\t", strip.white=TRUE)

mintotal = 1500
minexpr = 5
minnumber = 1
maxexpr = 500
dodownsample = FALSE
dsn = 1
rseed = 17000

for(file in annotations$fileprefix) {
  
  # Provided data has to be loaded
  setwd("D:/Userdata/jj.koning/MIKE/lymphnode-sc-transcriptomics-seperatescripts/data/1 - rawcounts")
  print(paste("Filtering File ", file, sep= ""))  
  data <-read.csv(paste(file,'.coutt.csv', sep=""), sep="\t", header=TRUE)
  
  # The first column contains unique gene ids and has to be assigned as rownames:
  rownames(data) <-data$GENEID
  
  # Remove gene IDs (of the non-endogenous spike-in RNAs), starting with ERCC: 
  data <-data[grep("ERCC",rownames(data),invert=TRUE),-1]
  
  # Filtering of expression data (previously through RaceID)
  if ( ! is.numeric(mintotal) ) stop( "mintotal has to be a positive number" ) else if ( mintotal <= 0 ) stop( "mintotal has to be a positive number" )
  if ( ! is.numeric(minexpr) ) stop( "minexpr has to be a non-negative number" ) else if ( minexpr < 0 ) stop( "minexpr has to be a non-negative number" )
  if ( ! is.numeric(minnumber) ) stop( "minnumber has to be a non-negative integer number" ) else if ( round(minnumber) != minnumber | minnumber < 0 ) stop( "minnumber has to be a non-negative integer number" )
  if ( ! ( is.numeric(dodownsample) | is.logical(dodownsample) ) ) stop( "dodownsample has to be logical (TRUE/FALSE)" )
  if ( ! is.numeric(dsn) ) stop( "dsn has to be a positive integer number" ) else if ( round(dsn) != dsn | dsn <= 0 ) stop( "dsn has to be a positive integer number" )
  
  downsample <- function(x,n,dsn){
    x <- round( x[,apply(x,2,sum,na.rm=TRUE) >= n], 0)
    nn <- min( apply(x,2,sum) )
    for ( j in 1:dsn ){
      z  <- data.frame(GENEID=rownames(x))
      rownames(z) <- rownames(x)
      initv <- rep(0,nrow(z))
      for ( i in 1:dim(x)[2] ){
        y <- aggregate(rep(1,nn),list(sample(rep(rownames(x),x[,i]),nn)),sum)
        na <- names(x)[i]
        names(y) <- c("GENEID",na)
        rownames(y) <- y$GENEID
        z[,na] <- initv
        k <- intersect(rownames(z),y$GENEID)
        z[k,na] <- y[k,na]
        z[is.na(z[,na]),na] <- 0
      }
      rownames(z) <- as.vector(z$GENEID)
      ds <- if ( j == 1 ) z[,-1] else ds + z[,-1]
    }
    ds <- ds/dsn  # Removed +0.1
    return(ds)
  }
  
  ndata <- data[,apply(data,2,sum,na.rm=TRUE) >= mintotal]
  if ( dodownsample){
      set.seed(rseed)
      ndata <- downsample(data,n=mintotal,dsn=dsn)
  }else{
      x <- ndata
      ndata <- as.data.frame( t(t(x)/apply(x,2,sum))*median(apply(x,2,sum,na.rm=TRUE))   ) # Removed +0.1
  }
  x <- ndata
  fdata <- x[apply(x>=minexpr,1,sum,na.rm=TRUE) >= minnumber,]
  x <- fdata
  fdata <- x[apply(x,1,max,na.rm=TRUE) < maxexpr,]
  
  # Write results
  fdata <- add_rownames(fdata, "GENEID")
  setwd("D:/Userdata/jj.koning/MIKE/lymphnode-sc-transcriptomics-seperatescripts/data/2 - filteredcounts")
  write_tsv(fdata, paste(file,'.coutt.csv', sep=""), append = FALSE)
  

}
