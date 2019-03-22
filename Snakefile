# Base scRNAseq analysis Snakemake that combines all other Snakefiles in ./workflows together. 

# Include common settings
include: "./workflows/common.snakefile"

# Modulization through Subworkflows: 
subworkflow star_index:
        workdir:
            "./workflows/star-index"
        snakefile:
            "./workflows/star-index/Snakefile"

subworkflow star_align:
        workdir:
            "./workflows/star-align"
        snakefile:
            "./workflows/star-align/Snakefile"

subworkflow map_barcodes:
        workdir:
            "./workflows/map-barcodes"
        snakefile:
            "./workflows/map-barcodes/Snakefile"


subworkflow quality_filters:
        workdir:
            "./workflows/quality-filters"
        snakefile:
            "./workflows/quality-filters/Snakefile"


subworkflow combine_plates:
        workdir:
            "./workflows/combine-plates"
        snakefile:
            "./workflows/combine-plates/Snakefile"


subworkflow raceid_stemid:
    workdir:
        "./workflows/raceid-stemid"
    snakefile:
        "./workflows/raceid-stemid/Snakefile"

