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

# Combine All Tables
tables <- lapply(annotations$fileprefix, function(x) {read_tsv(paste('1 - rawcounts/',x,'.coutt.csv',sep="")) %>% mutate('Plate'=x)}) %>% lapply(function(x) {gather(x, key='cellnr', value='count', -c('Plate','GENEID'))}) %>% bind_rows() %>% mutate("CellID" = paste(Plate,cellnr, sep = ".")) 
cells <- tables %>% select("Plate", "CellID") %>% distinct()
tables <- tables %>% select(-"cellnr", -"Plate") %>% spread(CellID,count)
tables[is.na(tables)] <- 0

# Write results
write.csv(tables, file = '3 - combinedcounts/DUITS__ALL.csv', row.names = FALSE)
write.csv(cells, file = '3 - combinedcounts/DUITS_CELLIDS.csv', row.names = FALSE)