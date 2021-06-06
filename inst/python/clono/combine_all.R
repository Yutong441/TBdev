#!/usr/bin/env Rscript
root_dir <- commandArgs(TRUE)[1]
condition <- commandArgs(TRUE)[2]
root <- paste (root_dir, condition, 'proc_data', sep='/')
#root <- '/mnt/c/Users/Yutong/Documents/bioinformatics/reproduction/results/imaging/data/clono/proc_data'

library (magrittr)
list.files (root, pattern='*_sum.csv', full.names=T) %>% as.list () %>%
        lapply (read.csv) -> dat_list

all_data <- do.call (rbind, dat_list)
all_data %>% dplyr::select (!dplyr::one_of(c('X.1', 'X', 'Unnamed..0'))) %>%
        write.csv (paste (root, 'all_sum.csv', sep='/'))
