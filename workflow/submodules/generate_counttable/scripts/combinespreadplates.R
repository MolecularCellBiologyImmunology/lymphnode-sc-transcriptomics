  ### Load / install required packages
.libPaths(snakemake@params$rpackagesfolders)
library(readr)
library(tidyr)
library(dplyr)

### Load Gene Annotation / Conversation Table
marttable <- read_tsv(paste(snakemake@input$conversiontable, "/genelist.tsv", sep = ""))

### Load Sample Annotation Table
annotations <- read_tsv(snakemake@input$annotations)
if (snakemake@config['istest'] == TRUE) {annotations <- filter(annotations, annotations$fileprefix %in% snakemake@params$samples)}

if (snakemake@params$samplescolumn != "~") {
### Combine & Spread Plates (Per Sample and All), Spread, Write Outputs
alltables <- data.frame()
variable = colnames(snakemake@params$samplescolumn)[1]
for(samplenumber in unique(snakemake@params$samplescolumn)) {
  print(paste("*** Combining Plates for Sample", samplenumber, sep=" "))
  annotationspersample <- filter(annotations, variable == samplenumber)
  sampleplates <- data.frame()
  for(prefix in annotationspersample$fileprefix) {
    print(paste("Converting File ", prefix, sep= ""))  
    data <- read_tsv(file.path(snakemake@input$inputfolder, prefix ,'counts.tsv'))
    data$geneid <- marttable[match(data$geneid, marttable$ensembl_gene_id),]$mgi_symbol
    data <- filter(data, !is.na(geneid))
    data <- filter(data, geneid != "")
    data <- mutate(data, "Plate" = prefix)
    data <- mutate(data, "CellID" = paste(prefix, cellbc, sep = "."))
    sampleplates <- rbind(sampleplates, data)
    alltables <- rbind(alltables, data)}
  print(paste("*** Writing Combined Plates for Sample", samplenumber,sep=" "))
  sampleplates <- dplyr::select(sampleplates, -"cellbc", -"Plate")
  sampleplates <- group_by(sampleplates, geneid, CellID)
  sampleplates <- summarise(sampleplates, reads = sum(reads))
  sampleplates <- spread(sampleplates, key = CellID, value = reads, fill = 0)
  write.csv(sampleplates, file = paste(snakemake@params$outputfolder, "/sample_", samplenumber, ".csv", sep=""), row.names = FALSE)}
} else {
  for(prefix in annotationspersample$fileprefix) {
    print(paste("Converting File ", prefix, sep= ""))  
    data <- read_tsv(file.path(snakemake@input$inputfolder, prefix ,'counts.tsv'))
    data$geneid <- marttable[match(data$geneid, marttable$ensembl_gene_id),]$mgi_symbol
    data <- filter(data, !is.na(geneid))
    data <- filter(data, geneid != "")
    data <- mutate(data, "Plate" = prefix)
    data <- mutate(data, "CellID" = paste(prefix, cellbc, sep = "."))
    alltables <- rbind(alltables, data)}
}

print("*** Writing CellID/Plate Lookup Table for All Cells")
allcells <- distinct(dplyr::select(alltables, "Plate", "CellID"))
write.csv(allcells, file = paste(snakemake@params$outputfolder,"cell_ids_plates_lookup_table.csv",sep="/"), row.names = FALSE)

print("*** Writing Count Table For all Plates Combined")
alltables <- dplyr::select(alltables, -"cellbc", -"Plate")
alltables <- group_by(alltables, geneid, CellID)
alltables <- summarise(alltables, reads = sum(reads))
alltables <- spread(alltables, key = CellID, value = reads, fill = 0)
write.csv(alltables, file = paste(snakemake@params$outputfolder,"ALL.csv",sep="/"), row.names = FALSE)

print("Finished.")
