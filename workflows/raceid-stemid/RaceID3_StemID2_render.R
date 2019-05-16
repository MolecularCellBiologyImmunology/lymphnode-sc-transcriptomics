# Automatically Determine Packages to download and install
if (!"rmarkdown" %in% installed.packages()) install.packages("rmarkdown", repos = "http://cran.us.r-project.org")

# Set Sample
sample = "ALL"

# Set working directory and intput/output locations
#workingdirectory = "C:/Users/Mike/Documents/WORK/Bioinformatics Project Internship/Scripts/seperate-scripts/lymphnode-sc-transcriptomics"
workingdirectory = "D:/Userdata/jj.koning/MIKE/Seperate Scripts/"
#workingdirectory = "D:/Documents/SCHOOL/VU/2017-2018 Master Year 2/Project/Seperate Scripts/lymphnode-sc-transcriptomics/"

scriptdirectory = paste(workingdirectory, "/scripts", sep="")
inputdata = paste(workingdirectory, "/data - douwe final/2 - combinedcounts/LNS_", sample, ".csv", sep="")
conversiontable = paste(workingdirectory, "/data - douwe final/genelist.tsv", sep = "")
outputdirectory = paste(workingdirectory, "/data - douwe final/3 - raceid3stemid2/", sample, sep="")



# Set parameters for Markdown Knitting
parameters = list(
  
  # Directories & Data
  scriptdirectory = scriptdirectory,
  inputdata = inputdata,
  conversiontable = conversiontable,
  outputdirectory = outputdirectory,
  
  #Filtering
  mintotal = 1500,
  minexpr = 5,
  minnumber = 1,
  LBatch = NULL,
  knn = FALSE,
  CGenes = c("Pcna","Mki67","Malat1","Hspa1a","Jun", "Fos", "Ptprc"),
  FGenes = NULL,
  ccor = 0.4,
  
  # RaceID
  maxclustnr = 50,
  bootnr = 50,
  
  # StemID
  RunStemID = FALSE,
  pdishuf=2000,
  scthr = 0.3
  )

# Run Initial Stem/RaceID Algorithm and Render Initial Report 
rmarkdown::render(paste(scriptdirectory, "RaceID3_StemID2_sample.Rmd", sep="/"), "html_document", output_dir = outputdirectory, params = parameters)

# Run Post-Analysis on specific Genes/Clusters based on Initial Report 
rmarkdown::render(paste(scriptdirectory, "RaceID3_StemID2_postanalysis.Rmd", sep="/"), "html_document", output_dir = outputdirectory)
