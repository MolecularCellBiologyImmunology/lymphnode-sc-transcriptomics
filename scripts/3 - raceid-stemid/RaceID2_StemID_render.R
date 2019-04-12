install.packages("rmarkdown")
require(rmarkdown)

# set working directory
setwd("C:/Users/Mike/Documents/WORK/Bioinformatics Project Internship/Scripts/seperate-scripts/lymphnode-sc-transcriptomics/scripts/3 - raceid-stemid")

# load class definition and functions
source("RaceID2_StemID_class.R")

# run stem/raceID and render report 
rmarkdown::render("RaceID2_StemID_sample.Rmd", "html_document", output_dir = "C:/Users/Mike/Documents/WORK/Bioinformatics Project Internship/Scripts/seperate-scripts/lymphnode-sc-transcriptomics/data/4 - raceidstemid")