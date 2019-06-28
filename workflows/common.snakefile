# Common Import Settings in the Workflows
import pandas as pd
from snakemake.utils import min_version
from snakemake.shell import shell
import pathlib
import datetime

##### Set minimum Snakemake Version #####

min_version("5.2.4")

##### Set Date and Time #####
x = datetime.datetime.now()
date_time = "-" + str(x.year) + "-" + str(x.month) + "-" + str(x.day) + "_" + str(x.hour) + "-" + str(x.minute) + "-" + str(x.second)

##### Load Config and Sample Sheets #####

configfile: config_file

##### Set Main Data Path #####

data = pathlib.Path(config['data'])

##### Annotations #####

# File with Sample Annotations
annotationfile = data / config['samples']
samples = pd.read_table(annotationfile).set_index('fileprefix', drop=False)
if config['istest']:
   samples = samples.loc[config['testset']]
samples.drop(columns='fileprefix')
fileprefixes = samples['fileprefix']

# File with Reference Genome Annotations
referencefile = data / config['references']
references = pd.read_table(referencefile).set_index('reference', drop=False)
reference = config['reference']

# Conversion file
conversionfolder = data / config['conversiontables'] / config['reference'] 

# File with Cell Barcodes
cellbcfile = data / config['celbc']

# Column defining samples
if config['samplescolumn'] != FALSE:
    samplescolumn = pd.read_table(annotationfile).set_index(config['samplescolumn'], drop=False)
    samplescolumn.drop(columns=config['samplescolumn'])
    samplescolumn = samplescolumn[config['samplescolumn']]
else:
    samplescolumn = FALSE

##### Other Paths #####

# Directory for Temporary Files, Intermediate Results
tmpstore = str(data / config['tmpstore'])

# Directory where output files will be stored
output = pathlib.Path(str(data / config['output']))
raceidoutputsbydate = str(data / config['output'] / 'raceid3stemid2') + date_time

# Directory Where the STAR Index Files for a Reference Genome are stored 
indexdir = pathlib.Path(str(data / config['index'] / config['reference']))

# Expected index files
indexfiles = ['SA','SAindex','chrLength.txt','chrName.txt','chrNameLength.txt']

# The Location of the GTF File of a Reference Genome
gtffile = pathlib.Path(str(data / config['refdir'] / config['reference'] / references.loc[config['reference'], 'gtffile']))

# The Location of the Fasta file of a Reference Genome
compressedfastafile = pathlib.Path(str(data / config['refdir'] / config['reference'] / references.loc[config['reference'], 'genomefile']))

# Locations of R packages
rpackagesfolders = config['rpackagesfolders']

##### Filter Settings #####

mintotal = config['mintotal']
minexpr = config['minexpr']
minnumber = config['minnumber']
LBatch = config['LBatch']
knn = config['knn']
CGenes = config['CGenes']
FGenes = config['FGenes']
ccor = config['ccor']

##### RaceID/StemID Settings #####

# RaceID3
maxclustnr = config['maxclustnr']
bootnr = config['bootnr']

# StemID2
RunStemID = config['RunStemID']
pdishuf= config['pdishuf']
scthr = config['scthr']
