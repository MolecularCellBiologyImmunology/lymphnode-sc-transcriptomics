# Automatically Determine Packages to download and install
list.of.packages <- c("rmarkdown")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, repos = "http://cran.us.r-project.org")


# Run Initial Stem/RaceID Algorithm and Render Initial Report 
rmarkdown::render("RaceID3_StemID2_sample.Rmd", "html_document", output_dir = snakemake@params[['outputdir']], params =  list(
  inputdata = snakemake@input[[]],
  outputdirectory,
  mintotal = snakemake@params[['mintotal']],
  minexpr = snakemake@params[['minexpr']],
  minnumber = snakemake@params[['minnumber']],
  maxexpr = snakemake@params[['maxexpr']],
  downsample = snakemake@params[['downsample']],
  sfn = snakemake@params[['sfn']],
  hkn = snakemake@params[['hkn']],
  dsn = snakemake@params[['dsn']],
  rseed = snakemake@params[['rseed']],
  CGenes = snakemake@params[['CGenes']],
  FGenes = snakemake@params[['FGenes']],
  ccor = snakemake@params[['ccor']],
  maxclustnr = snakemake@params[['maxclustnr']],
  bootnr = snakemake@params[['bootnr']],
  RunStemID = snakemake@params[['RunStemID']],
  pdishuf= snakemake@params[['pdishuf']],
  scthr = snakemake@params[['scthr']]))

# Run Post-Analysis on specific Genes/Clusters based on Initial Report 
rmarkdown::render("RaceID3_StemID2_postanalysis.Rmd", "html_document", output_dir = outputdirectory)
