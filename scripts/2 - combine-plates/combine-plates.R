# Install Packages (only first time)
install.packages("readr", "dplyr", "tidyr")

# Load required packages
library(readr)
library(tidyr)
library(dplyr)

# Paths
#setwd("C:/Users/Mike/Documents/WORK/Bioinformatics Project Internship/Scripts/seperate-scripts/lymphnode-sc-transcriptomics/data - douwe")
setwd("D:/Userdata/jj.koning/MIKE/Seperate Scripts/data - douwe")

# Annotations
annotations <- read.csv("annotations.tsv", sep = "\t")

# Combine All Tables
tables <- lapply(annotations$fileprefix, function(x) {read_tsv(paste('2 - spreadcounts/',x,'.tsv',sep="")) %>% mutate('Plate'=x)}) %>% lapply(function(x) {gather(x, key='cellnr', value='count', -c('Plate','geneid'))}) %>% bind_rows() %>% mutate("CellID" = paste(Plate,cellnr, sep = ".")) 
cells <- tables %>% select("Plate", "CellID") %>% distinct()
tables <- tables %>% select(-"cellnr", -"Plate") %>% spread(CellID,count)
tables[is.na(tables)] <- 0

# Write results
write.csv(tables, file = '3 - combinedcounts/LNS_ALL.csv', row.names = FALSE)
write.csv(cells, file = '3 - combinedcounts/CELLIDS_ALL.csv', row.names = FALSE)