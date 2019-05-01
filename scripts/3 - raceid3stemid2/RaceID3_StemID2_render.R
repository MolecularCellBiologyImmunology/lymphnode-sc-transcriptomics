list.of.packages <- c("rmarkdown")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, repos = "http://cran.us.r-project.org")

workingdirectory = "C:/Users/Mike/Documents/WORK/Bioinformatics Project Internship/Scripts/seperate-scripts/lymphnode-sc-transcriptomics"
#workingdirectory = "D:/Documents/SCHOOL/VU/2017-2018 Master Year 2/Project/Seperate Scripts/lymphnode-sc-transcriptomics/data")
#workingdirectory = "D:/Userdata/jj.koning/MIKE/Seperate Scripts/"

scriptdirectory = paste(workingdirectory, "scripts/3 - raceid3stemid2", sep="/")
outputdirectory = paste(workingdirectory, "data/4 - raceidstemid/W2", sep="/")
inputdata = paste(workingdirectory, 'data/3 - combinedcounts/LNS_W2_unfiltered_combined.csv', sep='/')

parameters = list(
  scriptdirectory = scriptdirectory,
  inputdata = inputdata,
  outputdirectory = outputdirectory,
  
  mintotal = 1500,
  minexpr = 5,
  minnumber = 1,
  maxexpr = 500,
  downsample = FALSE,
  sfn = FALSE,
  hkn = FALSE,
  dsn = 1,
  rseed = 17000,
  CGenes = c("Pcna","Mki67","Malat1","Hspa1a","Jun", "Fos"),
  FGenes = NULL,
  ccor = 0.4,
  
  maxclustnr = 30,
  bootnr = 50,
  
  RunStemID = TRUE,
  pdishuf=2000,
  scthr = 0.3)

# Run Initial Stem/RaceID Algorithm and Render Initial Report 
rmarkdown::render(paste(scriptdirectory, "RaceID3_StemID2_sample.Rmd", sep="/"), "html_document", output_dir = outputdirectory, params = parameters)

# Run Post-Analysis on specific Genes/Clusters based on Initial Report 
rmarkdown::render("scripts/3 - raceid3/RaceID3_StemID2_postanalysis.Rmd", "html_document", output_dir = "data/4 - raceidstemid/Duits")
