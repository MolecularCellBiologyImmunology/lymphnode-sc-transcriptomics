rmarkdown::render(
    'RaceID2_StemID_sample.Rmd', 
    input_file = snakemake@countscombined,
    output_file = snakemake@outputfile,
    output_dir = snakemake@outputdir,
)
