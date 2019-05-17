  ### Load / install required packages
#install.packages("readr", "dplyr", "tidyr", repos = "http://cran.us.r-project.org")
library(readr)
library(tidyr)
library(dplyr)

### Load Gene Annotation / Conversation Table
marttable <- read_tsv(paste(snakemake@input[[conversionfolder]], "/genelist.tsv", sep = ""))

### Load Sample Annotation Table
annotations <- read_tsv(snakemake@input[[annotations]])

### Combine & Spread Plates (Per Week and All), Spread, Write Outputs
alltables <- data.frame()
for(samplenumber in unique(annotations$time)) {
  print(paste("*** Combining Plates for Week", samplenumber, sep=" "))
  annotationspertime <- filter(annotations, time == samplenumber)
  sampleplates <- data.frame()
  for(prefix in annotationspertime$fileprefix) {
    print(paste("Converting File ", prefix, sep= ""))  
    data <- read_tsv(paste(snakemake@input[[inputfolder]], prefix ,'counts.tsv', sep="/"))
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
  write.csv(sampleplates, file = paste(snakemake@input[[outputfolder]], "/counttables_combined_sample_", samplenumber, ".csv", sep=""), row.names = FALSE)}

print("*** Writing CellID/Plate Lookup Table for All Cells")
allcells <- distinct(dplyr::select(alltables, "Plate", "CellID"))
write.csv(allcells, file = paste(snakemake@input[[outputfolder]],"cell_ids_plates_lookup_table.csv",sep="/"), row.names = FALSE)

print("*** Writing Count Table For all Plates Combined")
alltables <- dplyr::select(alltables, -"cellbc", -"Plate")
alltables <- group_by(alltables, geneid, CellID)
alltables <- summarise(alltables, reads = sum(reads))
alltables <- spread(alltables, key = CellID, value = reads, fill = 0)
write.csv(alltables, file = paste(snakemake@input[[outputfolder]],"counttables_combined_ALL.csv",sep="/"), row.names = FALSE)

print("Finished.")
