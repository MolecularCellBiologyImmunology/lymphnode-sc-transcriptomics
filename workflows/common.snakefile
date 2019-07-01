# Common Import Settings in the Workflows
import pandas as pd
from snakemake.utils import validate, min_version
from snakemake.shell import shell
import pathlib
import datetime

##### Set minimum Snakemake Version #####

min_version("5.2.4")

##### Set Date and Time #####
x = datetime.datetime.now()
date_time = "-" + str(x.year) + "-" + str(x.month) + "-" + str(x.day) + "_" + str(x.hour) + "-" + str(x.minute) + "-" + str(x.second)

##### Load Config and Sample Sheets #####
# We test whether it is defined, i.e. whether this is called from the 
# main workflow (then it's defined) or from a subworkflow (then it's 
# not defined) 
try:
    config_file
except NameError:
    config_file = None

if config_file is None:
    config_file = '../../snakemake_config.yaml'

configfile: config_file
validate(config, schema="schemas/config.schema.yaml")

##### Set Main Data Path #####

data = pathlib.Path(config['data'])

##### Annotations #####

# File with Sample Annotations
annotationfile = data / config['samples']
samples = pd.read_table(annotationfile).set_index('fileprefix', drop=False)
validate(samples, schema="schemas/annotations.schema.yaml")
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
if config['samplescolumn'] != '':
    samplescolumn = pd.read_table(annotationfile).set_index(config['samplescolumn'], drop=False)
    samplescolumn.drop(columns=config['samplescolumn'])
    samplescolumn = samplescolumn[config['samplescolumn']]
else:
    samplescolumn = ''

##### Other Paths #####

# Directory for Temporary Files, Intermediate Results
tmpstore = pathlib.Path(str(data / config['tmpstore']))

# Directory where output files will be stored
output = pathlib.Path(str(data / config['output']))
raceidoutputsbydate = pathlib.Path(str(data / config['output'] / 'raceid3stemid2') + date_time)

# Directory where the reference genome and gtf files are located
genomedir = data / config['refdir'] / config['reference']

# Directory where the STAR index files for a reference genome are located
config['indexdir'] = str(data / config['index'] / config['reference'])

# Paths of the reference files
# Fasta file
config['fastafile'] = str(genomedir /  references.loc[config['reference'], 'genome'])
# GTF file
config['gtffile'] = str(genomedir / references.loc[config['reference'], 'gtf'])

# Expected index files
# indexfiles = ['SA','SAindex','chrLength.txt','chrName.txt','chrNameLength.txt']

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
