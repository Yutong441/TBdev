#!/usr/bin/env Rscript
root <- commandArgs(TRUE)[1]
condition <- commandArgs(TRUE)[2]

#root <- '/mnt/c/Users/Yutong/Documents/bioinformatics/reproduction/results/imaging/'
data_dir <- paste (root, '/data/proc_data/', sep='')

tree <- read.csv (paste (data_dir, condition, '_.csv', sep=''))
DAPI <- read.csv (paste (data_dir, condition, '_nuc.csv', sep=''))
bright <- read.csv (paste (data_dir, condition, '_cyto.csv', sep=''))

devtools::load_all ('../../../', export_all=F)
plot_data <- preprocess (DAPI, bright, tree)
write.csv (plot_data, paste (data_dir, '/', condition, '_sum.csv',sep=''))
