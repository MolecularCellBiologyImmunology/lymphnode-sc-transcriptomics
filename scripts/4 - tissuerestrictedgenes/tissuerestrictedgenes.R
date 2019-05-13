install.packages("readr", " dplyr")
library(readr)
library(dplyr)

# Input Tissue Restricted Genes
#setwd("C:/Users/Mike/Documents/WORK/Bioinformatics Project Internship/Scripts/seperate-scripts/lymphnode-sc-transcriptomics/scripts/4 - tissuerestrictedgenes")
setwd("D:/Userdata/jj.koning/MIKE/Seperate Scripts/misc")
#setwd("D:/Documents/SCHOOL/VU/2017-2018 Master Year 2/Project/Seperate Scripts/lymphnode-sc-transcriptomics/scripts/4 - tissuerestrictedgenes")
y <- read_tsv("Sansom Supplemental_Table2.txt", skip=2)
y = y[-1,]
colnames(y)[colnames(y)=="Gene Symbol"] <- "mgi_symbol"
interesting = which(y$`GNF GeneAtlas Specificity` == "not present")
y = y[-interesting,]
interesting = which(y$`GNF GeneAtlas Specificity` == "not detected")
y = y[-interesting,]
interesting = which(y$`GNF GeneAtlas Specificity` == "NRE")
y = y[-interesting,]
output = data.frame()

# For Each cluster
#setwd("C:/Users/Mike/Documents/WORK/Bioinformatics Project Internship/Scripts/seperate-scripts/lymphnode-sc-transcriptomics/data/4 - raceidstemid/AllPlates")
setwd("D:/Userdata/jj.koning/MIKE/Seperate Scripts/data - douwe/4 - raceid3stemid2/ALL")
for (file in list.files(path = ".",pattern="cell_clust_diff_genes_")) {

  # Input RaceID StemID Output
  x <- read.csv(file, sep="\t", header = TRUE)
  
  # Filtered Join
  z <- inner_join(x,y, by = 'mgi_symbol')
  
  # Add Cluster Nr.
  file = sub(".xls","",sub("cell_clust_diff_genes_","",file))
  z <- mutate(z, "Cluster" = file)
  z = data.frame("Cluster" = z$Cluster, "Ensembl Code" = z$`Ensembl ID`, "MGI Symbol" = z$mgi_symbol, "Type Restricted" = x$`GNF GeneAtlas Specificity`, "Mean Outside Cluster" = z$mean.ncl, "Mean Inside Cluster" = z$mean.cl, "Factor Difference" = z$fc, "P Value" = z$pv, "Adjusted P Value" = z$padj, "Chrosome" = z$chromosome_name, "Annotation" = z$description)
  
  output <- rbind(output,z)
}

# Write Output
#setwd("C:/Users/Mike/Documents/WORK/Bioinformatics Project Internship/Scripts/seperate-scripts/lymphnode-sc-transcriptomics/data/5 - tissuerestrictedgenes/")
setwd("D:/Userdata/jj.koning/MIKE/Seperate Scripts/data - douwe/5 - tissuerestrictedgenes/")
output <- output[order(output["Cluster"], output["GNF GeneAtlas Specificity"], output["mgi_symbol"], decreasing = FALSE ),]
row.names(output) <- NULL 
write.table(output, 'tissuerestrictedgenes.tsv', sep = "\t", row.names = FALSE)
