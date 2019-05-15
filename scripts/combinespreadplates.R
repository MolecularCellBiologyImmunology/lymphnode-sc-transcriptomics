### Load / install required packages
#install.packages("readr", "dplyr", "tidyr", repos = "http://cran.us.r-project.org")
library(readr)
library(tidyr)
library(dplyr)

### Set Main Working Directory
setwd("C:/Users/Mike/Documents/WORK/Bioinformatics Project Internship/Scripts/seperate-scripts/lymphnode-sc-transcriptomics/data - douwe final")
#setwd("D:/Documents/SCHOOL/VU/2017-2018 Master Year 2/Project/Seperate Scripts/lymphnode-sc-transcriptomics/data - douwe final")
#setwd("D:/Userdata/jj.koning/MIKE/Seperate Scripts/data - douwe final")

### Load Gene Annotation / Conversation Table
marttable <- read_tsv("genelist.tsv")

### Load Sample Annotation Table
annotations <- read_tsv("annotations.tsv")

### Combine & Spread Plates (Per Week and All), Spread, Write Outputs
alltables <- data.frame()
for(samplenumber in unique(annotations$time)) {
  print(paste("*** Combining Plates for Week", samplenumber, sep=" "))
  annotationspertime <- filter(annotations, time == samplenumber)
  sampleplates <- data.frame()
  for(prefix in annotationspertime$fileprefix) {
    print(paste("Converting File ", prefix, sep= ""))  
    data <- read_tsv(paste("1 - rawcounts/", prefix ,'/counts.tsv', sep=""))
    data$geneid <- marttable[match(data$geneid ,marttable$ensembl_gene_id),]$mgi_symbol
    data <- filter(data, !is.na(geneid))
    data <- filter(data, geneid != "")
    data <- mutate(data, "Plate" = prefix)
    data <- mutate(data, "CellID" = paste(prefix, cellbc, sep = "."))
    sampleplates <- rbind(sampleplates, data)
    alltables <- rbind(alltables, data)}
  print(paste("*** Writing Combined Plates for Week", samplenumber,sep=" "))
  sampleplates <- dplyr::select(sampleplates, -"cellbc", -"Plate")
  sampleplates <- group_by(sampleplates, geneid, CellID)
  sampleplates <- summarise(sampleplates, reads = sum(reads))
  sampleplates <- spread(sampleplates, key = CellID, value = reads, fill = 0)
  write.csv(sampleplates, file = paste('2 - combinedcounts/LNS_W', samplenumber, ".csv", sep=""), row.names = FALSE)}

print("*** Writing CellID/Plate Lookup Table for All Cells")
allcells <- distinct(dplyr::select(alltables, "Plate", "CellID"))
write.csv(allcells, file = '2 - combinedcounts/CELLIDS_PLATES.csv', row.names = FALSE)

print("*** Writing Count Table For all Plates Combined")
alltables <- dplyr::select(alltables, -"cellbc", -"Plate")
alltables <- group_by(alltables, geneid, CellID)
alltables <- summarise(alltables, reads = sum(reads))
alltables <- spread(alltables, key = CellID, value = reads, fill = 0)
write.csv(alltables, file = '2 - combinedcounts/LNS_ALL.csv', row.names = FALSE)

print("Finished.")
