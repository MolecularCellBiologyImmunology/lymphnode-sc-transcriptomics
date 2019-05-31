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

subworkflow combine_spread_plates:
        workdir:
            "./workflows/combine-spread-plates"
        snakefile:
            "./workflows/combine-spread-plates/Snakefile"

subworkflow raceid_stemid:
    workdir:
        "./workflows/raceid-stemid"
    snakefile:
        "./workflows/raceid-stemid/Snakefile"


# Main workflow rule, determining output of workflow
rule all:
    input:
        index = star_index(expand("{indexdir}/{indexfile}", indexdir=indexdir, indexfile=indexfiles)),
        alignment = star_align(expand('{tmpstore}/star-align/{reference}/{fileprefix}/Aligned.out.bam', tmpstore=tmpstore, reference = reference, fileprefix=fileprefixes)),
        counts = map_barcodes(expand('{tmpstore}/mapping/{fileprefix}/counts.tsv', tmpstore=tmpstore, fileprefix=fileprefixes)),
        countscombined = combine_spread_plates(expand('{output}/combined_count_tables/counttables_combined_ALL.csv', output=output)),
        stemidreport = raceid_stemid(expand('{raceidoutputsbydate}/RaceID3_StemID2_sample.html', raceidoutputsbydate=raceidoutputsbydate)),
