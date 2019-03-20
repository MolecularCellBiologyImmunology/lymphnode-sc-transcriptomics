# Required Packages
install.packages("purrr")
install.packages("dplyr")
install.packages("tidyr")
library(readr)
library(purrr)
library(dplyr)
library(tidyr)


# Paths and Annotations         TODO: file.path()
setwd("C:/Users/Mike/Documents/WORK/Bioinformatics Project Internship/Scripts/lymphnode-sc-transcriptomics/workflows/combine-plates")
annotations <- read.csv("annotations.tsv", sep = "\t")
annotations <- filter(annotations, experiment == "timeseries", time == "2")
setwd("C:/Users/Mike/Documents/WORK/Bioinformatics Project Internship/C:/Users/Mike/Documents/WORK/Bioinformatics Project Internship/Data/research-sbi-mouse-lymph-node-stromal-cell-time-series\Heidelberg_count_tables_2018_Aug/Heidelberg_count_tables_2018_Aug")


### Old Way
# tables <- list()
# for (prefix in annotations$fileprefix) {
#   table <- as.data.frame(read.csv(paste(prefix,'.coutt.csv',sep=""), sep="\t", header=TRUE))
#   for (cell in colnames(table)) {
#     if (cell != "GENEID") colnames(table)[colnames(table)==cell] <- paste(prefix,cell,sep="_")
#   }
#   tables[[prefix]] <- table
# }
# tables <- reduce(tables, full_join, by="GENEID")


tables <- lapply(annotations$fileprefix, function(x) {read_tsv(paste(x,'.coutt.csv',sep="")) %>% mutate('Plate'=x)}) %>% lapply(function(x) {gather(x, key='cellnr', value='count', -c('Plate','GENEID'))}) %>% bind_rows() %>% mutate("CellID" = paste(Plate,cellnr, sep = ".")) 
cells <- tables %>% select("Plate", "CellID") %>% unique()
tables <- tables %>% select(-"cellnr", -"Plate") %>% spread(CellID,count)
tables[is.na(tables)] <- 0


# Write results
write.csv(tables, file = "C:/Users/Mike/Documents/WORK/Bioinformatics Project Internship/Scripts/lymphnode-sc-transcriptomics/workflows/combine-plates/normalized_counts.csv")
