# Input RaceID StemID Output

x <- read.csv("C:/Users/Mike/Documents/WORK/Bioinformatics Project Internship/Scripts/seperate-scripts/lymphnode-sc-transcriptomics/data/4 - raceidstemid/W2/cell_clust_diff_genes_cl_6",sep=",",header=TRUE)
rownames(x) <- x$GENEID

y <- read.csv("C:/Users/Mike/Documents/WORK/Bioinformatics Project Internship/Scripts/seperate-scripts/lymphnode-sc-transcriptomics/data/5 - tissuerestrictedgenes/Sansom Supplemental_Table2.xlsx",sep="\t",header=TRUE))
