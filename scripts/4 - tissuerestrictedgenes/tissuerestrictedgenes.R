install.packages("readr", " dplyr")
library(readr)
library(dplyr)

# Input Tissue Restricted Genes

setwd("C:/Users/Mike/Documents/WORK/Bioinformatics Project Internship/Scripts/seperate-scripts/lymphnode-sc-transcriptomics/scripts/4 - tissuerestrictedgenes")
y <- read_tsv("Sansom Supplemental_Table2.txt", skip=2)
y = y[-1,]
colnames(y)[colnames(y)=="Gene Symbol"] <- "GENEID"
interesting = which(y$`GNF GeneAtlas Specificity` == "not present")
y = y[-interesting,]
interesting = which(y$`GNF GeneAtlas Specificity` == "not detected")
y = y[-interesting,]
interesting = which(y$`GNF GeneAtlas Specificity` == "NRE")
y = y[-interesting,]
output = data.frame()

# For Each cluster

setwd("C:/Users/Mike/Documents/WORK/Bioinformatics Project Internship/Scripts/seperate-scripts/lymphnode-sc-transcriptomics/data/4 - raceidstemid/All Plates")

for (file in list.files(path = ".",pattern="cell_clust_diff_genes_cl_")) {

  # Input RaceID StemID Output
  x <- read.csv(file, sep="\t", header = TRUE)
  
  # Filtered Join
  z <- inner_join(x,y, by = 'GENEID')
  
  # Add Cluster Nr.
  file = substr(file,26,nchar(file)-4)
  z <- mutate(z, "Cluster" = file)
  z <- select(z, "Cluster", everything())
  
  output <- rbind(output,z)
}

setwd("C:/Users/Mike/Documents/WORK/Bioinformatics Project Internship/Scripts/seperate-scripts/lymphnode-sc-transcriptomics/data/5 - tissuerestrictedgenes/")
output <- output[order( output["Cluster"], output["GNF GeneAtlas Specificity"], output["GENEID"], decreasing = FALSE ),  ]
row.names(output) <- NULL 
write.table(output, 'tissuerestrictedgenes.txt', sep = "\t")
