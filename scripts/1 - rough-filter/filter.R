# Note: This Script is mostly based on the instruction manual of RaceID2 and StemID.
install.packages("dplyr")
install.packages("readr")
library(dplyr)
library(readr)

setwd("D:/Userdata/jj.koning/MIKE/seperatescripts v3/data")
#setwd("C:/Users/Mike/Documents/WORK/Bioinformatics Project Internship/Scripts/seperate-scripts/lymphnode-sc-transcriptomics/data")
annotations <- read.csv("annotations.tsv", sep = "\t", strip.white=TRUE)

mintotal = 1500
minexpr = 5
minnumber = 1
maxexpr = 500


for(file in annotations$fileprefix) {
  
  # Provided data has to be loaded
  setwd("D:/Userdata/jj.koning/MIKE/seperatescripts v3/data/1 - rawcounts")
  #setwd("C:/Users/Mike/Documents/WORK/Bioinformatics Project Internship/Scripts/seperate-scripts/lymphnode-sc-transcriptomics/data/1 - rawcounts")
  print(paste("Filtering File ", file, sep= ""))  
  data <-read.csv(paste(file,'.coutt.csv', sep=""), sep="\t", header=TRUE)
  
  # The first column contains unique gene ids and has to be assigned as rownames:
  rownames(data) <-data$GENEID
  
  # Remove gene IDs (of the non-endogenous spike-in RNAs) starting with ERCC: 
  data <-data[grep("ERCC",rownames(data),invert=TRUE),-1]
  
  # Rough filtering of expression data (to reduce table combining time)
  if ( ! is.numeric(mintotal) ) stop( "mintotal has to be a positive number" ) else if ( mintotal <= 0 ) stop( "mintotal has to be a positive number" )
  if ( ! is.numeric(minexpr) ) stop( "minexpr has to be a non-negative number" ) else if ( minexpr < 0 ) stop( "minexpr has to be a non-negative number" )
  if ( ! is.numeric(minnumber) ) stop( "minnumber has to be a non-negative integer number" ) else if ( round(minnumber) != minnumber | minnumber < 0 ) stop( "minnumber has to be a non-negative integer number" )

  
  data <- data[,apply(data,2,sum,na.rm=TRUE) >= mintotal]
  data <- data[apply(data>=minexpr,1,sum,na.rm=TRUE) >= minnumber,]
  fdata <- data[apply(data,1,max,na.rm=TRUE) < maxexpr,]
  
  # Write results
  fdata <- add_rownames(fdata, "GENEID")
  setwd("D:/Userdata/jj.koning/MIKE/seperatescripts v3/data/2 - filteredcounts")
  #setwd("C:/Users/Mike/Documents/WORK/Bioinformatics Project Internship/Scripts/seperate-scripts/lymphnode-sc-transcriptomics/data/2 - filteredcounts")
  write_tsv(fdata, paste(file,'.coutt.csv', sep=""), append = FALSE)


}
