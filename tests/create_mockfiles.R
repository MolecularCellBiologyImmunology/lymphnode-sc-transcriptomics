# creating mock sequence files for testing purposes
library(readr)
an <- read_tsv('/home/douwe/data/mouse-lymph-node-stromal-cell-time-series/annotations_2018_Aug/annotations.csv')

mocklocation <- '/home/douwe/data/mouse-lymph-node-stromal-cell-time-series/sequences_2018_Aug'
mockfiles <- file.path(mocklocation, paste0(an$filePrefix, '_R1.fastq.gz'))
lapply(mockfiles, function(x){
  cmd <- paste('touch', x)
  system(cmd)
})
