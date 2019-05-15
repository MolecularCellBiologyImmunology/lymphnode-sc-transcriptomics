# Automatically Determine Packages to download and install
if (!"rmarkdown" %in% installed.packages()) install.packages("rmarkdown", repos = "http://cran.us.r-project.org")


# Set working directory and intput/output locations
#workingdirectory = "C:/Users/Mike/Documents/WORK/Bioinformatics Project Internship/Scripts/seperate-scripts/lymphnode-sc-transcriptomics"
#workingdirectory = "D:/Userdata/jj.koning/MIKE/Seperate Scripts/"
workingdirectory = "D:/Documents/SCHOOL/VU/2017-2018 Master Year 2/Project/Seperate Scripts/lymphnode-sc-transcriptomics/"
scriptdirectory = paste(workingdirectory, "scripts/3 - raceid3stemid2", sep="/")
inputdata = paste(workingdirectory, 'data - douwe final/2 - combinedcounts/LNS_ALL.csv', sep='/')
outputdirectory = paste(workingdirectory, "data - douwe final/3 - raceid3stemid2/ALL", sep="/")



# Set parameters for Markdown Knitting
parameters = list(
  
  # Directories & Data
  scriptdirectory = scriptdirectory,
  inputdata = inputdata,
  outputdirectory = outputdirectory,
  
  #Filtering
  mintotal = 1500,
  minexpr = 1,
  minnumber = 3,
  LBatch = NULL,
  knn = FALSE,
  CGenes = c("Pcna","Mki67","Malat1","Hspa1a","Jun", "Fos", "Ptprc"),
  FGenes = NULL,
  ccor = 0.4,
  
  # RaceID
  maxclustnr = 30,
  bootnr = 50,
  
  # StemID
  RunStemID = TRUE,
  pdishuf=2000,
  scthr = 0.3
  )

# Run Initial Stem/RaceID Algorithm and Render Initial Report 
rmarkdown::render(paste(scriptdirectory, "RaceID3_StemID2_sample.Rmd", sep="/"), "html_document", output_dir = outputdirectory, params = parameters)

# Run Post-Analysis on specific Genes/Clusters based on Initial Report 
rmarkdown::render(paste(scriptdirectory, "RaceID3_StemID2_postanalysis.Rmd", sep="/"), "html_document", output_dir = outputdirectory)
