# ==========embed the data in GPLVM==========
root_dir <- '/mnt/c/Users/Yutong/Documents/bioinformatics/reproduction/'
root <- paste (root_dir, 'results/', sep='/')
merge_dir <- paste (root, 'XLYBPZ_Dylan_dir', sep='/')
save_dir <- paste (root, 'manuscript/figure2', sep='/')
all_data <- readRDS(paste (merge_dir, 'revision_all.rds', sep='/') )

gp_dir <- paste (save_dir, 'GPLVM_embedding', sep='/')
if (!dir.exists (gp_dir)) {dir.create (gp_dir)}

# load python files
reticulate::use_condaenv("ptime")
#script_root <- system.file("python/", package='TBdev')
script_root <- paste (root_dir, 'TBdev/inst/python/', sep='')
reticulate::source_python (paste (script_root, 'GPLVM/GrandPrix.py', sep='/'))
library (TBdev)

# preprocess data
tb_data <- all_data [, ! (all_data$revised %in% c(CT$in_vitro_cells, CT$non_emb_lineage,
                   'Oocyte', 'Zy', '2C', '4C', '8C', 'PE', 'EPI', 'PSA-EPI' ))]
tb_data <- tb_data [, tb_data$date != 'in_vitro']
exp_mat <- save_to_csv (tb_data, select_genes=4000, directory=save_dir)
Y <- scale(t(exp_mat))

# run GPLVM
reticulate::py_set_seed (100)
Grand_out = fit_model (data=Y, # rows are cells and columns are genes
                       n_latent_dims=3L, #here 3 latent dimensions are computed
                       n_inducing_points=10L, #number of inducing points
                       latent_prior_var=1., 
                       maxitr=100L, #maximum epoch of iteration
                       disp=T,  return_pred=T,
                       kernel=list ('name'='RBF', 'ls'=1.0, 'var'=1.0)
                       # Use RBF kernel
)

tb_data <- RunGPLVM (tb_data, Grand_out[[1]][[1]])
plot_dim_red (tb_data, group.by = c('broad_type', 'date'), DR='gplvm',
              dims=c(1,2,3), all_theta=50)

saveRDS (Grand_out, paste (save_dir, 'GPLVM_embedding/GPLVM_embedding.rds', sep='/'))
GPLVM2csv (Grand_out, Y, paste (save_dir, 'GPLVM_embedding', sep='/'))

