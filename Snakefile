# Base scRNAseq analysis Snakemake that combines all other Snakefiles in ./workflows together. 
import easygui

##### Load Config and Sample Sheets #####
print("Please choose a configuration file to run Snakemake with.")
config_file = easygui.fileopenbox()

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

# subworkflow filters:
#         workdir:
#             "./workflows/quality-filters"
#         snakefile:
#             "./workflows/quality-filters/Snakefile"

subworkflow combineplates:
        workdir:
            "./workflows/combine-plates"
        snakefile:
            "./workflows/combine-plates/Snakefile"

subworkflow raceid_stemid:
    workdir:
        "./workflows/raceid-stemid"
    snakefile:
        "./workflows/raceid-stemid/Snakefile"


# Main workflow rule, determining output
rule all:
    input:
        index = star_index(expand("{indexdir}/{indexfile}", indexdir=indexdir, indexfile=indexfiles)),
        alignment = star_align(expand('{tmpstore}/star-align/{fileprefix}/Aligned.out.bam', tmpstore=tmpstore, fileprefix=fileprefixes)),
        countsfile = map_barcodes(expand('{tmpstore}/mapping/{fileprefix}/counts.tsv',  tmpstore=tmpstore, fileprefix=fileprefixes)),
        countsreport = map_barcodes(expand('{tmpstore}/mapping/{fileprefix}/report.html', tmpstore=tmpstore, fileprefix=fileprefixes)),
        # filtered = filters(expand('{tmpstore}/quality_filters/{fileprefix}/counts_filtered.tsv', tmpstore=tmpstore, fileprefix=fileprefixes)),
        combined = combineplates(expand('{output}/counts_combined.tsv', output=output)),
        stemidreport = raceid_stemid(expand('{output}/RaceID2_StemID_results.html', output=output)),
        clustercells = raceid_stemid(expand('{output}/cell_clust.xlx', output=output)),
        clustergenes = raceid_stemid(expand('{output}/cell_clust_diff_genes_cl_*', output=output))
