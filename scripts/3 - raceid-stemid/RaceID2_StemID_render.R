install.packages("rmarkdown")
require(rmarkdown)

# set working directory
setwd("D:/Documents/SCHOOL/VU/2017-2018 Master Year 2/Project/Seperate Scripts/lymphnode-sc-transcriptomics/scripts/3 - raceid-stemid")

# load class definition and functions
source("RaceID2_StemID_class.R")

# Run Initial Stem/RaceID Algorithm and Render Initial Report 
rmarkdown::render("RaceID2_StemID_sample.Rmd", "html_document", output_dir = "D:/Documents/SCHOOL/VU/2017-2018 Master Year 2/Project/Seperate Scripts/lymphnode-sc-transcriptomics/data/4 - raceidstemid/AllPlates")

# Run Post-Analysis on specific Genes/Clusters based on Initial Report 
rmarkdown::render("RaceID2_StemID_postanalysis.Rmd", "html_document", output_dir = "D:/Documents/SCHOOL/VU/2017-2018 Master Year 2/Project/Seperate Scripts/lymphnode-sc-transcriptomics/data/4 - raceidstemid/AllPlates")