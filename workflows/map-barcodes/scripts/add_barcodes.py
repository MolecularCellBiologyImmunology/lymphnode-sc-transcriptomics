# snakemake workflow version
# add barcodes (umi and cell) to a sam file

# input[fastqfile]: the fastq (fastq.gz) file with the barcodes
# input[samfile]: the sam file
# output[0]: the sam file with the barcodes added as fields at the end of the lines 

import sys
from gatb import Bank
import sqlite3
from common.utils import chunks

r = []
print("Reading source file into in-memory database of barcodes: {}".format(snakemake.input['fastqfile']))
con = sqlite3.connect(":memory:")
cur = con.cursor()
print("")
cur.execute('CREATE TABLE barcodes (seqid TEXT, celbc TEXT, umi TEXT)')
# reading from file and writing to database in chunks to save memory
fastq_parser = Bank(snakemake.input.fastqfile[0])
# TODO: can this be normalized to snakemake.input.fastqfile since the path was 
# changed from a pathlib object to a string object??
for chunk in chunks(fastq_parser, 10000):
    r = []
    for seq in chunk:
        sequence = seq.sequence.decode("utf-8")
        umi = sequence[0:6]
        cel = sequence[6:12]
        seqid = seq.comment.decode("utf-8").split(" ")[0]
        r.append((seqid, cel, umi))
    cur.executemany('INSERT INTO barcodes VALUES (?,?,?)', r)
r = None
print("Creating index on read indentifiers")
cur.execute('CREATE UNIQUE INDEX seqidx ON barcodes (seqid)')

print("Writing output file: {}".format(snakemake.output[0]))
# TODO: can this be normalized to snakemake.output since the path was 
# changed from a pathlib object to a string object??
fi = open(snakemake.input.samfile, 'r') 
# TODO: can this be normalized to snakemake.output since the path was 
# changed from a pathlib object to a string object??
fo = open(snakemake.output[0], 'w', buffering=32768)
for line in fi:
    if line.startswith('@'):
        fo.write(line)
    else:
        line = line[:-1] # removing '\n' at end of line
        elements = line.split('\t')
        cur.execute('SELECT celbc, umi FROM barcodes WHERE seqid=:seqid', {'seqid': elements[0]})
        bcs = cur.fetchone()
        elements.append('bc:Z:' + bcs[0])
        elements.append('um:Z:' + bcs[1])
        fo.write('\t'.join(elements) + '\n')
fi.close()
fo.close()
