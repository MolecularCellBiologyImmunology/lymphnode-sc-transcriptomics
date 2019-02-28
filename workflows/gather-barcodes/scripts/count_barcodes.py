# snakemake workflow version
# count barcodes (umi and cell) from a fastq sequence file

# input[0]: the fastq (fastq.gz) file
# input[1]: the file with known cell barcodes
# output[0]: the file with all code combinations counted
# output[1]: the file with only known cell codes

import sys
from gatb import Bank
import csv

# read the known cell codes
celbc = []
print("Reading celBC file: {}".format(snakemake.input[1]))
with open(snakemake.input[1]) as tsvfile:
    tsvreader = csv.reader(tsvfile, delimiter='\t')
    for row in tsvreader:
        celbc = celbc + [row[1]]

r = {}
print("Reading source file: {}".format(snakemake.input[0]))
fi = open(snakemake.output[1], 'w')
fastq_parser = Bank(snakemake.input[0])
for seq in fastq_parser:
    sequence = seq.sequence.decode("utf-8")
    umi = sequence[0:6]
    cel = sequence[6:12]
    seqid = seq.comment.decode("utf-8").split(" ")[0]
    if cel in r.keys():
        if umi in r[cel].keys():
            r[cel][umi] += 1
        else:
            r[cel][umi] = 1
    else:
        r[cel] = {umi:1}
    if cel in celbc:
        fi.write("{id},{cel},{umi}\n".format(id=seqid, cel=cel, umi=umi))
fi.close()

print("Opening summary file: {}".format(snakemake.output[0]))
with open(snakemake.output[0], 'w') as fo:
    for cel in r.keys():
        for umi in r[cel].keys():
            fo.write("{fileprefix},{cel},{umi},{count}\n".format(
                fileprefix=snakemake.wildcards['fileprefix'], 
                cel=cel, 
                umi=umi, 
                count=str(r[cel][umi])))
