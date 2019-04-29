# Common Import Settings in the Workflows
import pandas as pd
from snakemake.utils import min_version
from snakemake.shell import shell
import pathlib

##### Set minimum Snakemake Version #####

min_version("5.2.4")

##### Load Config and Sample Sheets #####

configfile: config_file

##### Set Data Path #####

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

# File with Cell Barcodes
cellbcfile = data / config['celbc']

##### Other Paths #####

# Directory for Temporary Files, Intermediate Results
tmpstore = str(data / config['tmpstore'])

# Directory where output files will be stored
output = str(data / config['output'])

# Directory Where the STAR Index Files for a Reference Genome are stored 
indexdir = str(data / config['index'] / config['reference'])

# Expected index files
indexfiles = ['SA','SAindex','chrLength.txt','chrName.txt','chrNameLength.txt']

# The Location of the GTF File of a Reference Genome
gtffile = str(data / config['refdir'] / config['reference'] / references.loc[config['reference'], 'gtffile'])

# The Location of the Fasta file of a Reference Genome
fastafile = str(data / config['refdir'] / config['reference'] / references.loc[config['reference'], 'genomefile'])

##### Filter Settings #####

mintotal = config['mintotal']
minexpr = config['minexpr']
minnumber = config['minnumber']
maxexpr = config['maxexpr']
dodownsample = config['dodownsample']
dsn = config['dsn']
rseed = config['rseed']

##### RaceID/StemID Settings #####
