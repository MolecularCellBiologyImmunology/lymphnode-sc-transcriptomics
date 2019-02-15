# Note: This Script is mostly based on the instruction manual of RaceID2 and StemID.
# First this script, then RaceID2_StemID_sample.R should be run.

# Install & load required packages
# install.packages(c("kulife", "XML"))
library(kulife)
library(XML)

########## TODO: Import output file of generate_counts here ##########
# Provided data has to be loaded
data <-read.csv("testfile.xls", sep="\t", header=TRUE)

# The first column contains unique gene ids and has to be assigned as rownames:
rownames(data) <-data$GENEID

# Remove gene IDs (of the non-endogenous spike-in RNAs), starting with “ERCC”: 
data <-data[grep("ERCC",rownames(data),invert=TRUE),-1]

# Filtering of expression data (previously through RaceID)
# TODO: Copy the main changes this function makes from class script
# data <- filterdata(sc, mintotal=3000, minexpr=5, minnumber=1, maxexpr=500, downsample=TRUE, dsn=1, rseed=17000)

# Write results
write.csv(data, file = "testfile_cleaned.csv", row.names=FALSE,col.names=TRUE,sep="\t",quote=FALSE)
