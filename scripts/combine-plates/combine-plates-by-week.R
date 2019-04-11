# Install Packages (only first time)
install.packages("readr")
install.packages("dplyr")
install.packages("tidyr")

# Load required packages
library(readr)
library(tidyr)
library(dplyr)

# Get file prefixes for main loop
setwd("D:/Userdata/jj.koning/MIKE/lymphnode-sc-transcriptomics-seperatescripts/data")
annotations <- read.csv("annotations.tsv", sep = "\t")

# Combine plates for each experiment time point
for(filter in unique(annotations$time)) {

  print(paste("Combining Plates for Week",filter,sep=" "))
  
  # Paths and Annotations
  setwd("D:/Userdata/jj.koning/MIKE/lymphnode-sc-transcriptomics-seperatescripts/data")
  annotations <- read.csv("annotations.tsv", sep = "\t")
  annotations <- filter(annotations, time == filter)
  setwd("D:/Userdata/jj.koning/MIKE/lymphnode-sc-transcriptomics-seperatescripts/data/1 - rawcounts")

  # Combining Tables per Timeframe
  tables <- lapply(annotations$fileprefix, function(x) {read_csv(paste(x,'.coutt.csv',sep="")) %>% mutate('Plate'=x)}) %>% lapply(function(x) {gather(x, key='cellnr', value='count', -c('Plate','GENEID'))}) %>% bind_rows() %>% mutate("CellID" = paste(Plate,cellnr, sep = ".")) 
  cells <- tables %>% select("Plate", "CellID") %>% unique()
  tables <- tables %>% select(-"cellnr", -"Plate") %>% spread(CellID,count)
  tables[is.na(tables)] <- 0
  
  # Write results
  setwd("D:/Userdata/jj.koning/MIKE/lymphnode-sc-transcriptomics-seperatescripts/data/2 - combinedcounts")
  write.csv(tables, file = paste('LNS_W',filter,'.coutt.csv',sep=""), row.names = FALSE)
  write.csv(tables, file = paste('LNS_W',filter,'.cellplateids.csv',sep=""), row.names = FALSE)

}