.libPaths(c("C:/Users/Mike/Documents/R/win-library/3.6", "C:/Program Files/R/R-3.6.0/library"))
library(rmarkdown)  

# Set Sample
sample = "ALL"

# Set parameters for Markdown Knitting
parameters = list(
  
  # Directories & Data
  inputdata       = snakemake@input$countscombined,
  conversiontable = paste(snakemake@input$conversionfolder, "/genelist.tsv", sep = ""),
  outputdirectory = snakemake@params$outputdirectory,
  
  #Filtering
  mintotal    = snakemake@params$mintotal,
  minexpr     = snakemake@params$minexpr,
  minnumber   = snakemake@params$minnumber,
  LBatch      = snakemake@params$LBatch,
  knn         = snakemake@params$knn,
  CGenes      = snakemake@params$CGenes,
  FGenes      = snakemake@params$FGenes,
  ccor        = snakemake@params$ccor,
  
  # RaceID
  maxclustnr  = snakemake@params$maxclustnr,
  bootnr      = snakemake@params$bootnr,
  
  # StemID
  RunStemID   = snakemake@params$RunStemID,
  pdishuf     = snakemake@params$pdishuf,
  scthr       = snakemake@params$scthr
  )

# Run Initial Stem/RaceID Algorithm and Render Initial Report 
rmarkdown::render("scripts/RaceID3_StemID2_sample.Rmd", "html_document", output_dir = snakemake@params$outputdirectory, params = parameters)

# Run Post-Analysis on specific Genes/Clusters based on Initial Report 
#rmarkdown::render(paste("RaceID3_StemID2_postanalysis.Rmd"), "html_document", output_dir = snakemake@params[[outputdirectory]])