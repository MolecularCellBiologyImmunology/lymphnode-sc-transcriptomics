install.packages("rmarkdown")
require(rmarkdown)

# set working directory
setwd("D:/Documents/SCHOOL/VU/2017-2018 Master Year 2/Project/Seperate Scripts/lymphnode-sc-transcriptomics/scripts/3 - raceid-stemid")

# load class definition and functions
source("RaceID2_StemID_class.R")

# run stem/raceID and render report 
rmarkdown::render("RaceID2_StemID_sample.Rmd", "html_document", output_dir = "D:/Documents/SCHOOL/VU/2017-2018 Master Year 2/Project/Seperate Scripts/lymphnode-sc-transcriptomics/data/4 - raceidstemid")