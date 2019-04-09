# Required Packages
install.packages("purrr")
install.packages("dplyr")
install.packages("tidyr")
install.packages("readr")
library(readr)
library(purrr)
library(tidyr)
library(dplyr)


for(filter in unique(annotations$time)) {

  # Paths and Annotations
  setwd("D:/Documents/SCHOOL/VU/2017-2018 Master Year 2/Project/Seperate Scripts/lymphnode-sc-transcriptomics/data")
  annotations <- read.csv("annotations.tsv", sep = "\t")
  annotations <- filter(annotations, experiment == "timeseries", time == filter)
  
  setwd("D:/Documents/SCHOOL/VU/2017-2018 Master Year 2/Project/Seperate Scripts/lymphnode-sc-transcriptomics/data/2 - filteredcounts")
  tables <- lapply(annotations$fileprefix, function(x) {read_tsv(paste(x,'.coutt.csv',sep="")) %>% mutate('Plate'=x)}) %>% lapply(function(x) {gather(x, key='cellnr', value='count', -c('Plate','GENEID'))}) %>% bind_rows() %>% mutate("CellID" = paste(Plate,cellnr, sep = ".")) 
  cells <- tables %>% select("Plate", "CellID") %>% unique()
  tables <- tables %>% select(-"cellnr", -"Plate") %>% spread(CellID,count)
  tables[is.na(tables)] <- 0
  
  # Write results
  setwd("D:/Documents/SCHOOL/VU/2017-2018 Master Year 2/Project/Seperate Scripts/lymphnode-sc-transcriptomics/data/3 - combinedcounts")
  write.csv(tables, file = paste('LNS_W',filter,'.coutt.csv',sep=""))

}
