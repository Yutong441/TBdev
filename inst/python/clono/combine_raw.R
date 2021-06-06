#!/usr/bin/env Rscript
root <- commandArgs(TRUE)[1]
condition <- commandArgs(TRUE)[2]

data_dir <- paste (root, '/proc_data/', sep='')
DAPI <- read.csv (paste (data_dir, condition, '_.csv', sep=''), stringsAsFactors=F)
DAPI$condition <- basename (DAPI$condition)
sub1 <- gsub('_nuc_seg.tiff$', '', DAPI$series)
DAPI$series <- gsub ('^.*_', '', sub1)
write.csv (DAPI, paste (data_dir, '/', condition, '_sum.csv',sep=''))
