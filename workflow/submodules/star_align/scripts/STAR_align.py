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

fileprefix = os.path.splitext(os.path.basename(input_str))[0]
outprefix = os.path.join(os.path.dirname(outputdir) + "/" + fileprefix + "/")

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