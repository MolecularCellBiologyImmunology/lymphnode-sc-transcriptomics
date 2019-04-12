# Install Packages (only first time)
install.packages("readr")
install.packages("dplyr")
install.packages("tidyr")

# Load required packages
library(readr)
library(tidyr)
library(dplyr)

# Get file prefixes for main loop
setwd("D:/Userdata/jj.koning/MIKE/seperatescripts v3/data")
annotations <- read.csv("annotations.tsv", sep = "\t")

# Combine plates for each experiment time point
for(filter in unique(annotations$time)) {

  print(paste("Combining Plates for Week",filter,sep=" "))
  
  # Paths and Annotations      
  setwd("D:/Userdata/jj.koning/MIKE/seperatescripts v3/data")
  annotations <- read.csv("annotations.tsv", sep = "\t")
  annotations <- filter(annotations, time == filter)
  setwd("D:/Userdata/jj.koning/MIKE/seperatescripts v3/data/2 - filteredcounts")
  
  # Combine All Tables
  tables <- lapply(annotations$fileprefix, function(x) {read_tsv(paste(x,'.coutt.csv',sep="")) %>% mutate('Plate'=x)}) %>% lapply(function(x) {gather(x, key='cellnr', value='count', -c('Plate','GENEID'))}) %>% bind_rows() %>% mutate("CellID" = paste(Plate,cellnr, sep = ".")) 
  cells <- tables %>% select("Plate", "CellID") %>% unique()
  tables <- tables %>% select(-"cellnr", -"Plate") %>% spread(CellID,count)
  tables[is.na(tables)] <- 0
  
  # Write results
  setwd("D:/Userdata/jj.koning/MIKE/seperatescripts v3/data//3 - combinedcounts")
  write.csv(tables, file = paste('LNS_W',filter,'_combined.csv',sep=""), row.names = FALSE)
  write.csv(cells, file = paste('LNS_W',filter,'_cellplateids.csv',sep=""), row.names = FALSE)

}