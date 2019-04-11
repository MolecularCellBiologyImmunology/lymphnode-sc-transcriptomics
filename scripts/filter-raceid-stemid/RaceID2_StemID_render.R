install.packages("rmarkdown")
require(rmarkdown)

# set working directory
setwd("D:/Userdata/jj.koning/MIKE/lymphnode-sc-transcriptomics-seperatescripts/scripts/filter-raceid-stemid")

# load class definition and functions
source("RaceID2_StemID_class.R")

# run stem/raceID and render report 
rmarkdown::render("RaceID2_StemID_sample.R", "pdf_document")