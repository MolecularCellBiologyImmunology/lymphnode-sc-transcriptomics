# generate simple table of unique umi per true cell barcode and gene
# NOTE: this script is not functional. Please see the R-script with the 
# same name for an alternative

import sqlite3
import csv
from common.utils import chunks

con = sqlite3.connect(":memory:")
cur = con.cursor()
r = []
print("Loading cell barcode file to in-memory database")
cur.execute('DROP TABLE IF EXISTS cellbc')
cur.execute('CREATE TABLE cellbc (bcnr INTEGER, cellbc TEXT PRIMARY KEY)')
cb = open(snakemake.input['cellbcfile'], 'r')
for line in cb:
    r.append(line.strip().split('\t'))
cb.close()
cur.executemany('INSERT INTO cellbc VALUES (?, ?)', r)

print("Loading hq mapped feature table to in-memory database: {}".format(snakemake.input['featuretable']))
cur.execute('DROP TABLE IF EXISTS featuremaps')
cur.execute('CREATE TABLE featuremaps (cellbc TEXT, umi TEXT, geneid TEXT)')
# reading from file and writing to database in chunks to save memory
with open(snakemake.input['featuretable'], 'r') as ft:
    dictreader = csv.DictReader(ft, dialect='excel-tab')
    for chunk in chunks(dictreader, 10000):
        r = []
        for rw in chunk:
            r.append((rw['cellbc'], rw['umi'], rw['geneid']))
        cur.executemany('INSERT INTO featuremaps VALUES (?,?,?)', r)
r = None

print("Creating indexes on features")
cur.execute('CREATE INDEX cellbcidx ON featuremaps (cellbc)')
cur.execute('CREATE INDEX umiidx ON featuremaps (umi)')
cur.execute('CREATE INDEX geneididx ON featuremaps (geneid)')

cur.executescript("""
SELECT cellbc, geneid, count(umi) FROM cellbc 
LEFT JOIN featuremaps USING (cellbc) 
GROUP BY cellbc, geneid
LIMIT 100
""")
print(cur.fetchall())
cur.close()
con.close()

# print("Joining tables and counting unique cellbc-umi-geneid combinations")

# cur.executescript("""

# """)
