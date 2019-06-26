# Base scRNAseq analysis Snakemake that combines all other Snakefiles in ./workflows together. 

##### Load Config and Sample Sheets #####
import easygui
print("Please choose a configuration file to run Snakemake with.")
config_file = easygui.fileopenbox()

# Include common settings
include: "./workflows/common.snakefile"

# Modulization through Subworkflows: 
subworkflow star_index:
        snakefile:
            "./workflows/star-index/Snakefile"

subworkflow star_align:
        snakefile:
            "./workflows/star-align/Snakefile"

subworkflow map_barcodes:
        snakefile:
            "./workflows/map-barcodes/Snakefile"

subworkflow combine_spread_plates:
        snakefile:
            "./workflows/combine-spread-plates/Snakefile"

subworkflow raceid_stemid:
    snakefile:
        "./workflows/raceid-stemid/Snakefile"

# Main workflow rule, determining output of workflow
rule all:
    input:
        module1 = star_index(expand("{indexdir}/{indexfile}", indexdir=indexdir, indexfile=indexfiles)),
        module2 = star_align(expand('{tmpstore}/star-align/{reference}/{fileprefix}/Aligned.out.bam', tmpstore=tmpstore, reference = reference, fileprefix=fileprefixes)),
        module3 = map_barcodes(expand('{tmpstore}/mapping/{fileprefix}/counts.tsv', tmpstore=tmpstore, fileprefix=fileprefixes)),
        module4 = combine_spread_plates(expand('{output}/combined_count_tables/counttables_combined_ALL.csv', output=output)),
        module5 = raceid_stemid(expand('{raceidoutputsbydate}/RaceID3_StemID2_sample.html', raceidoutputsbydate=raceidoutputsbydate)),
