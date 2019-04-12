# Install Packages (only first time)
install.packages("readr")
install.packages("dplyr")
install.packages("tidyr")

# Load required packages
library(readr)
library(tidyr)
library(dplyr)

# Paths and Annotations      
setwd("C:/Users/Mike/Documents/WORK/Bioinformatics Project Internship/Scripts/seperate-scripts/lymphnode-sc-transcriptomics/data")
annotations <- read.csv("annotations.tsv", sep = "\t")
setwd("C:/Users/Mike/Documents/WORK/Bioinformatics Project Internship/Scripts/seperate-scripts/lymphnode-sc-transcriptomics/data/1 - rawcounts")

# Combine All Tables
tables <- lapply(annotations$fileprefix, function(x) {read_tsv(paste(x,'.coutt.csv',sep="")) %>% mutate('Plate'=x)}) %>% lapply(function(x) {gather(x, key='cellnr', value='count', -c('Plate','GENEID'))}) %>% bind_rows() %>% mutate("CellID" = paste(Plate,cellnr, sep = ".")) 
cells <- tables %>% select("Plate", "CellID") %>% unique()
tables <- tables %>% select(-"cellnr", -"Plate") %>% spread(CellID,count)
tables[is.na(tables)] <- 0

# Write results
  setwd("C:/Users/Mike/Documents/WORK/Bioinformatics Project Internship/Scripts/seperate-scripts/lymphnode-sc-transcriptomics/data/3 - combinedcounts")
  write.csv(tables, file = paste('normalized_counts.csv',sep=""), row.names = FALSE)
  write.csv(cells, file = paste('cellplateids.csv',sep=""), row.names = FALSE)