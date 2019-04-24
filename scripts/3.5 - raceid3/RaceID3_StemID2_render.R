list.of.packages <- c("rmarkdown")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, repos = "http://cran.us.r-project.org")



setwd("D:/Documents/SCHOOL/VU/2017-2018 Master Year 2/Project/Seperate Scripts/lymphnode-sc-transcriptomics/data")
annotations <- read.csv("annotations.tsv", sep = "\t")
#setwd("D:/Documents/SCHOOL/VU/2017-2018 Master Year 2/Project/Seperate Scripts/lymphnode-sc-transcriptomics/scripts/3.5 - raceid3")
setwd("D:/Documents/SCHOOL/VU/2017-2018 Master Year 2/Project/Seperate Scripts/lymphnode-sc-transcriptomics/scripts/3.5 - raceid3")



# Run Initial Stem/RaceID Algorithm and Render Initial Report 
rmarkdown::render("RaceID3_StemID2_sample.Rmd", "html_document", output_dir = "D:/Documents/SCHOOL/VU/2017-2018 Master Year 2/Project/Seperate Scripts/lymphnode-sc-transcriptomics/data/4 - raceidstemid/AllPlates")

# Run Post-Analysis on specific Genes/Clusters based on Initial Report 
rmarkdown::render("RaceID3_StemID2_postanalysis.Rmd", "html_document", output_dir = "D:/Documents/SCHOOL/VU/2017-2018 Master Year 2/Project/Seperate Scripts/lymphnode-sc-transcriptomics/data/4 - raceidstemid/AllPlates")
