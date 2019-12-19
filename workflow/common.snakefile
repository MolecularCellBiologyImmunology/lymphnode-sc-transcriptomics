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
date_time = "-" + str(x.year) + "-" + str(x.month) + "-" + str(x.day) + \
    "_" + str(x.hour) + "-" + str(x.minute) + "-" + str(x.second)

##### Load Config and Sample Sheets #####
# We test whether it is defined, i.e. whether this is called from the 
# main workflow (then it's defined) or from a subworkflow (then it's 
# not defined) 
try:
    config_file
except NameError:
    config_file = None

if config_file is None:
    config_file = '../snakemake_config.yaml'

configfile: config_file
validate(config, schema="schemas/config.schema.yaml")

##### Set Main Data Path #####

datadir = pathlib.Path(config['datadir'])

##### Annotations #####

# File with Sample Annotations
annotationfile = str(datadir / config['repo'] / config['samples'])
samples = pd.read_table(annotationfile).set_index('fileprefix', drop=False)
validate(samples, schema="schemas/annotations.schema.yaml")
if config['istest']:
   samples = samples.loc[config['testset']]
samples.drop(columns='fileprefix')
fileprefixes = samples['fileprefix']

# File with Reference Genome Annotations
referencefile = str(datadir / config['repo'] / config['references'])
references = pd.read_table(referencefile).set_index('reference', drop=False)
validate(references, schema='schemas/references.schema.yaml')

# Conversion file
conversionfolder = datadir / config['conversiontables'] / config['reference'] 

# File with Cell Barcodes
config['cellbcfile'] = str(datadir / config['repo'] / config['celbc'])

# Constructing directory with sequences
config['seqdir'] = str(datadir / config['repo'] / config['sequences']) 

# Column defining samples
if config['samplescolumn'] != '':
    samplescolumn = pd.read_table(annotationfile).set_index(config['samplescolumn'], drop=False)
    samplescolumn.drop(columns=config['samplescolumn'])
    samplescolumn = samplescolumn[config['samplescolumn']]
else:
    samplescolumn = ''

##### Other Paths #####

# Directory for Temporary Files, Intermediate Results
tmppath = str(datadir / config['repo'] / config['tmpstore'])

# Directory where output files will be stored
outputpath = pathlib.Path(str(datadir / config['output']))
raceidoutputsbydate = pathlib.Path(str(datadir / config['output'] / 'raceid3stemid2') + date_time)

# Constructing main directory where the reference genomes and gtf files are
# located
config['genomedir'] = str(datadir / config['refdir']) # / config['reference'] ?

# Constructing main directory where the STAR index files are located
config['indexdir'] = str(datadir / config['indexdir']) # / config['reference'] ?

# Names of the reference files
# Fasta file
config['fastafile'] = references.loc[config['reference'], 'genome']
# GTF file
config['gtffile'] = references.loc[config['reference'], 'gtf']

# Expected index files
indexfiles = ['SA','SAindex','chrLength.txt','chrName.txt','chrNameLength.txt']

# Locations of R packages
rpackagesfolders = config['rpackagesfolders']
