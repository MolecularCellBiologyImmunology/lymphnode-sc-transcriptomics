# Install Packages (only first time)
install.packages("readr")
install.packages("dplyr")
install.packages("tidyr")

# Load required packages
library(readr)
library(tidyr)
library(dplyr)

# Paths and Annotations      
setwd("D:/Userdata/jj.koning/MIKE/Seperate Scripts/data")
annotations <- read.csv("annotations.tsv", sep = "\t")
setwd("D:/Userdata/jj.koning/MIKE/Seperate Scripts/data/1 - rawcounts")

# Combine All Tables
tables <- lapply(annotations$fileprefix, function(x) {read_tsv(paste(x,'.coutt.csv',sep="")) %>% mutate('Plate'=x)}) %>% lapply(function(x) {gather(x, key='cellnr', value='count', -c('Plate','GENEID'))}) %>% bind_rows() %>% mutate("CellID" = paste(Plate,cellnr, sep = ".")) 
cells <- tables %>% select("Plate", "CellID") %>% distinct()
tables <- tables %>% select(-"cellnr", -"Plate") %>% spread(CellID,count)
tables[is.na(tables)] <- 0

# Write results
setwd("D:/Userdata/jj.koning/MIKE/Seperate Scripts/data/3 - combinedcounts")
write.csv(tables, file = paste('LNS_W_ALL.csv',sep=""), row.names = FALSE)
write.csv(cells, file = paste('LNS_W_ALL_CELLIDS.csv',sep=""), row.names = FALSE)