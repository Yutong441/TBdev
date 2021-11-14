# ==========Obtain STREAM embedding==========
# need to run 1GPLVM.R first
root_dir <- '/mnt/c/Users/Yutong/Documents/bioinformatics/reproduction/'
root <- paste (root_dir, 'results/', sep='/')
save_dir <- paste (root, 'manuscript/figure2/', sep='/')
stream_dir <- paste (save_dir, 'STREAM', sep='/')
if (!dir.exists (stream_dir)) {dir.create (stream_dir)}

reticulate::use_condaenv("scrna")
script_root <- system.file("python/", package='TBdev')
#script_root <- paste (root_dir, 'TBdev/inst/python/', sep='')
reticulate::source_python (paste (script_root, 'STREAM/utils.py', sep='/'))

reticulate::py_set_seed (100)
show_trajectory (save_dir)
# branch1: stem, branch2: EVT, branch3: STB
branch_dict <- list (c('S0', 'S1'), c('S2', 'S1'), c('S3', 'S1'))
names (branch_dict) <- paste ('branch', 1:3, sep='')

reticulate::py_set_seed (100)
save_trajectory (save_dir, stem='S0', branch_dict=branch_dict)
