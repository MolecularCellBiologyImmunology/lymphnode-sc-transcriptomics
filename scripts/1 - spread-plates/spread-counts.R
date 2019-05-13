# Install Packages (only first time)
install.packages("readr", "dplyr", "tidyr")
source("https://bioconductor.org/biocLite.R")
biocLite(pkgs=c("scran","DESeq2","biomaRt"), suppressUpdates = TRUE)

# Load required packages
library(readr)
library(tidyr)
library(dplyr)
library(biomaRt)

# Paths
setwd("C:/Users/Mike/Documents/WORK/Bioinformatics Project Internship/Scripts/seperate-scripts/lymphnode-sc-transcriptomics/data - douwe")
#setwd("D:/Userdata/jj.koning/MIKE/Seperate Scripts/data - douwe")

# Create Mart Table
mart <- useMart(biomart = "ensembl", dataset = "mmusculus_gene_ensembl")
marttable <- biomaRt::getBM(attributes = c("ensembl_gene_id", "mgi_symbol", "chromosome_name", "description"), mart = mart)
write_tsv(marttable, "genelist.tsv")


# Annotations
annotations <- read_tsv("annotations.tsv")

# Script
for(prefix in annotations$fileprefix) {
  print(paste("Gathering File ", prefix, sep= ""))  
  data <- read_tsv(paste("1 - rawcounts/", prefix ,'/counts.tsv', sep="")) 
  data <- filter(data, !is.na(geneid))
  data <- spread(data, key = cellbc, value = reads, fill = "0")
  genes <- data$geneid 
  data$geneid <- marttable[match(genes,marttable$ensembl_gene_id),]$mgi_symbol
  data <- write_tsv(data, paste("2 - spreadcounts/", prefix, '.tsv', sep=""))
}
