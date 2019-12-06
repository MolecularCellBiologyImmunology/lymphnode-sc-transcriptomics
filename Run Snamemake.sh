echo "Starting scRNA-seq analysis snakemake workflow."
cd "./scRNA Workflow/"
CALL conda.bat activate scRNA
Snakemake
pause
