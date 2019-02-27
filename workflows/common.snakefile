# Common settings in the workflows
import pandas as pd
from snakemake.utils import min_version
from snakemake.shell import shell
import pathlib

##### set minimum snakemake version #####
min_version("5.2.4")

##### load config and sample sheets #####

configfile: "../../config.yaml"

##### data path #####

data = pathlib.Path(config['data'])

##### annotations #####

# file with sample annotations
annotationfile = data / config['samples']
samples = pd.read_table(annotationfile).set_index('fileprefix', drop=False)
if config['istest']:
   samples = samples.loc[config['testset']]

# file with reference genome annotations
referencefile = data / config['references']
references = pd.read_table(referencefile).set_index('reference', drop=False)

# file with cell barcodes
cellbcfile = data / config['celbc']

##### other paths #####

# directory for temporary files, intermediate results
tmpstore = str(data / config['tmpstore'])
# directory where the STAR index files for a reference genome are stored 
indexdir = str(data / config['index'] / config['reference'])
# the location of the gtf file of a reference genome
gtffile = str(data / config['refdir'] / config['reference'] / references.loc[config['reference'], 'gtffile'])
# the location of the fasta file of a reference genome
# TODO: repair this
fastafile = str(data / config['refdir'] / config['reference'] / references.loc[config['reference'], 'genomefile'])
