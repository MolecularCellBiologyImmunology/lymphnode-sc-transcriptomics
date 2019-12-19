rmarkdown::render(
    'scripts/generate_counts.Rmd', 
    output_file = snakemake@params[['outputfile']],
    output_dir = snakemake@params[['outputdir']],
    params = list(
        featuretable = snakemake@input[['featuretable']],
        cellbcfile = snakemake@input[['cellbcfile']],
        countsfile = snakemake@output[['countsfile']],
        samplename = snakemake@wildcards[['fileprefix']],
        UMIcorrection = snakemake@params[['UMIcorrection']],
        UMIcardinality = snakemake@params[['UMIcardinality']]
    )
)