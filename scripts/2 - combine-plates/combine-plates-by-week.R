# Install Packages (only first time)
install.packages("readr")
install.packages("dplyr")
install.packages("tidyr")

# Load required packages
library(readr)
library(tidyr)
library(dplyr)

# Get file prefixes for main loop
setwd("D:/Documents/SCHOOL/VU/2017-2018 Master Year 2/Project/Seperate Scripts/lymphnode-sc-transcriptomics/data")
annotations <- read.csv("annotations.tsv", sep = "\t")

# Combine plates for each experiment time point
for(filter in unique(annotations$time)) {

  print(paste("Combining Plates for Week",filter,sep=" "))
  
  # Paths and Annotations      
  annotations <- read.csv("annotations.tsv", sep = "\t")
  annotations <- filter(annotations, time == filter)

  # Combine All Tables
  tables <- lapply(annotations$fileprefix, function(x) {read_tsv(paste('1 - rawcounts/',x,'.coutt.csv',sep="")) %>% mutate('Plate'=x)}) %>% lapply(function(x) {gather(x, key='cellnr', value='count', -c('Plate','GENEID'))}) %>% bind_rows() %>% mutate("CellID" = paste(Plate,cellnr, sep = ".")) 
  cells <- tables %>% select("Plate", "CellID") %>% distinct()
  tables <- tables %>% select(-"cellnr", -"Plate") %>% spread(CellID,count)
  tables[is.na(tables)] <- 0
  
  # Write results
  write.csv(tables, file = paste('3 - combinedcounts/LNS_W',filter,'_unfiltered_combined.csv',sep=""), row.names = FALSE)
  #write.csv(cells, file = paste('3 - combinedcounts/LNS_W',filter,'_cellplateids.csv',sep=""), row.names = FALSE)

}