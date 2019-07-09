### Set Main Working Directory
setwd("C:/Users/Mike/Documents/WORK/Bioinformatics Project Internship/Scripts/seperate-scripts/lymphnode-sc-transcriptomics/data - douwe final")
#setwd("D:/Documents/SCHOOL/VU/2017-2018 Master Year 2/Project/Seperate Scripts/lymphnode-sc-transcriptomics/data - douwe final")
#setwd("D:/Userdata/jj.koning/MIKE/Seperate Scripts/data - douwe final")

source("https://bioconductor.org/biocLite.R")
#biocLite(pkgs=c("scran","DESeq2","biomaRt", "readr"), suppressUpdates = TRUE)
library(biomaRt)
library(readr)

mart <- useMart(biomart = "ensembl", dataset = "mmusculus_gene_ensembl")
marttable <- biomaRt::getBM(attributes = c("ensembl_gene_id", "mgi_symbol", "chromosome_name", "description"), mart = mart)
write_tsv(marttable, "genelist.tsv")