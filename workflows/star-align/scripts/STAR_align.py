# adapted from https://snakemake-wrappers.readthedocs.io/en/stable/wrappers/star/align.html
# for single-end sequences

# __author__ = "Johannes Köster"
# __copyright__ = "Copyright 2016, Johannes Köster"
# __email__ = "koester@jimmy.harvard.edu"
# __license__ = "MIT"

import os
from snakemake.shell import shell
# from snakemake.utils import format

extra = snakemake.params.get("extra", "")
log = snakemake.log_fmt_shell(stdout=True, stderr=True)

fq = snakemake.input.get("fq")
assert fq is not None, "input-> fq is a required input parameter"
input_str = fq

# if fq.endswith(".gz"):
#     readcmd = "--readFilesCommand unpigz -c"
# else:
#     readcmd = ""

outprefix = os.path.dirname(snakemake.output[0]) + "/"

# print("Will execute the following statement:\n")
# print(format(
#     "STAR "
#     "{extra} "
#     "--runThreadN {snakemake.threads} "
#     "--genomeDir {snakemake.params.indexdir} "
#     "--readFilesIn {input_str} "
#     "--readFilesCommand unpigz -c "
#     "--outSAMtype BAM Unsorted "
#     "--outFileNamePrefix {outprefix} "
#     "--outStd Log "
#     "{log}"))

shell(
    "STAR "
    "{extra} "
    "--runThreadN {snakemake.threads} "
    "--genomeDir {snakemake.params.indexdir} "
    "--readFilesIn {input_str} "
    "--readFilesCommand unpigz -c "
    "--outSAMtype BAM Unsorted "
    "--outFileNamePrefix {outprefix} "
    "--outStd Log "
    "{log}")