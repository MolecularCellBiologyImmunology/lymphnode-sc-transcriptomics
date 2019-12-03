library(readr)
library(dplyr)

mapfile <- '/home/douwe/data/temporary/mapping/mouse-grcm38/LNS_W0_P2_7/HqFeatureCounts.tsv'
celbcfile <- '/home/douwe/data/mouse-lymph-node-stromal-cell-time-series/annotations_common/celseq_barcodes.192.txt'

d <- read_tsv(mapfile, comment='@', na='NA')
bc <- read_tsv(celbcfile, col_names = c('bcnr','cellbc'))

valid_celbc_counts <- d %>%
  group_by(cellbc) %>%
  summarise(count = n()) %>%
  left_join(bc) %>%
  group_by(is.na(bcnr)) %>%
  summarise(total=sum(count))

counts_per_cell_and_gene <- bc %>%
  left_join(d) %>%
  distinct(cellbc, geneid, umi) %>%
  group_by(cellbc, geneid) %>%
  summarise(reads = n())

reads_per_cell <- bc %>%
  left_join(d) %>%
  distinct(cellbc, geneid, umi) %>%
  group_by(cellbc) %>%
  summarise(reads = n())

genes_per_umi <- d %>%
  filter(mapquality=='255' & isprimarymap=='True') %>%
  distinct(celbc,umi,geneid) %>%
  group_by(celbc, umi) %>%
  summarise(genes = n())

replicates_per_umi <- d %>%
  filter(mapquality=='255' & isprimarymap=='True') %>%
  group_by(celbc, geneid, umi) %>%
  summarise(replicates = n()) %>%
  arrange(celbc, geneid, umi)

genes_per_cel <- d %>%
  filter(mapquality=='255' & isprimarymap=='True') %>%
  distinct(celbc, geneid) %>%
  group_by(celbc) %>%
  summarise(genespercel = n())
