install.packages("rmarkdown")
require(rmarkdown)

# set working directory
setwd("D:/Userdata/jj.koning/MIKE/seperatescripts v3/scripts/3 - filter-raceid-stemid")

# load class definition and functions
source("RaceID2_StemID_class.R")

# run stem/raceID and render report 
rmarkdown::render("RaceID2_StemID_sample.Rmd", "html_document", output_dir = "D:/Userdata/jj.koning/MIKE/seperatescripts v3/data/4 - raceidstemid")