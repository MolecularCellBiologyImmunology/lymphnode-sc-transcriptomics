list.of.packages <- c("rmarkdown")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, repos = "http://cran.us.r-project.org")

#setwd("D:/Documents/SCHOOL/VU/2017-2018 Master Year 2/Project/Seperate Scripts/lymphnode-sc-transcriptomics/data")
setwd("D:/Userdata/jj.koning/MIKE/Seperate Scripts/")

# Run Initial Stem/RaceID Algorithm and Render Initial Report 
rmarkdown::render("scripts/3 - raceid3/RaceID3_StemID2_sample.Rmd", "html_document", output_dir = "data/4 - raceidstemid/Duits")

# Run Post-Analysis on specific Genes/Clusters based on Initial Report 
rmarkdown::render("scripts/3 - raceid3/RaceID3_StemID2_postanalysis.Rmd", "html_document", output_dir = "data/4 - raceidstemid/Duits")
