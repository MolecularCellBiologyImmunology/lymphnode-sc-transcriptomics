# Install Packages (only first time)
install.packages("readr", "dplyr", "tidyr")

# Load required packages
library(readr)
library(tidyr)
library(dplyr)

# Get file prefixes for main loop
setwd("C:/Users/Mike/Documents/WORK/Bioinformatics Project Internship/Scripts/seperate-scripts/lymphnode-sc-transcriptomics/data - douwe")
#setwd("D:/Documents/SCHOOL/VU/2017-2018 Master Year 2/Project/Seperate Scripts/lymphnode-sc-transcriptomics/data - douwe")
#setwd("D:/Userdata/jj.koning/MIKE/Seperate Scripts/data - douwe")

# Combine plates for each experiment time point
for(filter in unique(annotations$time)) {

  print(paste("Combining Plates for Week",filter,sep=" "))
  
  # Paths and Annotations      
  annotations <- read.csv("annotations.tsv", sep = "\t")
  annotations <- filter(annotations, time == filter)

  # Combine All Tables
  tables <- lapply(annotations$fileprefix, function(x) {read_tsv(paste('2 - spreadcounts/',x,'.tsv',sep="")) %>% mutate('Plate'=x)}) %>% lapply(function(x) {gather(x, key='cellnr', value='count', -c('Plate','geneid'))}) %>% bind_rows() %>% mutate("CellID" = paste(Plate,cellnr, sep = ".")) 
  cells <- tables %>% select("Plate", "CellID") %>% distinct()
  tables <- tables %>% select(-"cellnr", -"Plate") %>% spread(CellID,count)
  tables[is.na(tables)] <- 0
  
  # Write results
  write.csv(tables, file = paste('3 - combinedcounts/LNS_W',filter,'_combined.csv',sep=""), row.names = FALSE)
}