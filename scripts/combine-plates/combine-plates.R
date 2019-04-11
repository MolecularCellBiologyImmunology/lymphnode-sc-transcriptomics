# Required Packages
#install.packages("readr")
#install.packages("dplyr")
#install.packages("tidyr")
library(readr)
library(tidyr)
library(dplyr)

# Paths and Annotations      
setwd("D:/Userdata/jj.koning/MIKE/lymphnode-sc-transcriptomics-seperatescripts/data")
annotations <- read.csv("annotations.tsv", sep = "\t")



setwd("D:/Userdata/jj.koning/MIKE/lymphnode-sc-transcriptomics-seperatescripts/data/2 - filteredcounts")


tables <- lapply(annotations$fileprefix, function(x) {read_tsv(paste(x,'.coutt.csv',sep="")) %>% mutate('Plate'=x)}) %>% lapply(function(x) {gather(x, key='cellnr', value='count', -c('Plate','GENEID'))}) %>% bind_rows() %>% mutate("CellID" = paste(Plate,cellnr, sep = ".")) 
cells <- tables %>% select("Plate", "CellID") %>% unique()
tables <- tables %>% select(-"cellnr", -"Plate") %>% spread(CellID,count)
tables[is.na(tables)] <- 0


# Write results
write.csv(tables, file = "D:/Userdata/jj.koning/MIKE/lymphnode-sc-transcriptomics-seperatescripts/data/3 - combinedcounts/normalized_counts.csv", row.names = FALSE)
