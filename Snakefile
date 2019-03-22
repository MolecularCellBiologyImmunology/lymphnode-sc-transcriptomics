# Base scRNAseq analysis Snakemake that combines all other Snakefiles in ./workflows together. 

# Include common settings
include: "./workflows/common.snakefile"

# Modulization through Subworkflows: 
subworkflow star-index
    workdir:
        workdir:
            "./workflows/star-index"
        snakefile:
            "./workflows/star-index/Snakefile"

subworkflow star-align
    workdir:
        workdir:
            "./workflows/star-align"
        snakefile:
            "./workflows/star-align/Snakefile"

subworkflow map-barcodes
    workdir:
        workdir:
            "./workflows/map-barcodes"
        snakefile:
            "./workflows/map-barcodes/Snakefile"


subworkflow quality-filters
    workdir:
        workdir:
            "./workflows/quality-filters"
        snakefile:
            "./workflows/quality-filters/Snakefile"


subworkflow combine-plates
    workdir:
        workdir:
            "./workflows/combine-plates"
        snakefile:
            "./workflows/combine-plates/Snakefile"


subworkflow raceid-stemid
    workdir:
        workdir:
            "./workflows/raceid-stemid"
        snakefile:
            "./workflows/raceid-stemid/Snakefile"

