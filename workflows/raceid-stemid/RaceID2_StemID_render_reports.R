rmarkdown::render(
    'RaceID2_StemID_sample.Rmd', 
    output_file = snakemake@params[['outputfile']],
    output_dir = snakemake@params[['outputdir']],
    params = list(

    )
)
