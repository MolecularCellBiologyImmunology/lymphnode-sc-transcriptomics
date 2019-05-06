# Install Packages (only first time)
install.packages("readr", "dplyr", "tidyr")

# Load required packages
library(readr)
library(tidyr)
library(dplyr)

# Paths
setwd("C:/Users/Mike/Documents/WORK/Bioinformatics Project Internship/Scripts/seperate-scripts/lymphnode-sc-transcriptomics/data - douwe")
#setwd("D:/Userdata/jj.koning/MIKE/Seperate Scripts/data - douwe")

# Annotations
annotations <- read_tsv("annotations.tsv", sep = "\t", strip.white=TRUE)

# Script
for(prefix in annotations$fileprefix) {
  print(paste("Gathering File ", prefix, sep= ""))  
  data <- read_tsv(paste("1 - rawcounts/", prefix ,'/counts.tsv', sep="")) 
  data <- spread(data, key = cellbc, value = reads, fill = "0") 
  data <- write_tsv(data, paste("2 - spreadcounts/", prefix, '.tsv', sep=""))
}
