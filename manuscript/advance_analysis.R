# This script demonstrates how to run WGCNA, quadratic programming, cell-cell
# pairwise correlation, pathway, cell cycle analysis and diffusion pseudotime.

root_dir <- '/mnt/c/Users/Yutong/Documents/bioinformatics/reproduction/'
setwd (paste (root_dir, 'SingleCellR/utils', sep='/'))
library (tidyverse)
library (Seurat)

ML <- modules::use ('marker_list.R')
DEG <- modules::use ('DE_gene.R')
KG <- modules::use ('KEGG_path.R')

# load data
root <- paste (root_dir, 'results/', sep='/')
merge_dir <- paste (root, 'XLYBPZ_Dylan_dir', sep='/')
save_dir <- paste (root, 'manuscript/figure2', sep='/')
sup_save_dir <- paste (root, 'manuscript/figureS2', sep='/')
x <- load (paste (merge_dir, 'final_merged.Robj', sep='/') )
all_data <- get (x)

# ----------pathway analysis----------
# obtain the module scores for each pathway for each cell
# Again due to the computational intensive nature of this calculation, we
# make it compulsory to save the object. 
# You may input a datarame in `all_path` containing the columns 'pathway' and
# 'kegg_id'. Otherwise, we have compiled a list of pathways in `KG$kg_pathway`
path_seurat <- KG$get_module_score (all_data, 
                                      all_path=KG$kg_pathway, 
                                      save_path=paste (sup_save_dir, 'Data_module_scores.csv', sep='/'),
                                      append_meta=T #return a Seurat object with metadata
)

# visualise the pathway scores in heatmap
DEG$seurat_heat (path_seurat, 
                 group.by=c('broad_type','date'),
                 normalize=T
)

# alternatively, visualise it with violin plot
DEG$seurat_violin (path_seurat, features=rownames (path_seurat), 
                   group.by='broad_type', 
                   slot_data='counts', #the data is stored in this slot
                   num_col=6 
)

# ----------cell cycle analysis----------
CC <- modules::use ('cycle_analysis.R')
exp_mat_cc <- as.matrix (GetAssayData (all_data, assay='RNA', slot='data'))
# calculate the phase score according to Macosco (2015)
ans <- CC$get_phase_score (exp_mat_cc)

# visualise using heatmap
CC$make_cycle_heat (ans, all_data, feature='broad_type')

# correlate cell cycle with PC
CC$phase_PC_plot (ans, all_data, corr_feature='broad_type', 
                  num_dim=1, # correlation with the PC1
)

# ----------WGCNA----------
# This function implements most of the functionality of WGCNA. This includes
# selection of an appropriate power for scale-free network, using this weight
# power to calculate connectivity matrix, and hierarchical clustering of the
# matrix. The eigen-gene of each cluster in each cell is returned.
# Owing to the various results it generates, we make it compulsory to store the
# results in a directory.
WG <- modules::use ('WGCNA_utils.R')
datME <- WG$find_eigengene (all_data, sup_save_dir, 
                            cluster_num='all' #use the entire dataset
)

# read the stored cluster information
color_row <- read.csv ( paste ( sup_save_dir, 'WGCNA/module_genes.csv' , sep='/'), row.names=1)
gene_list <- lapply (as.list (colnames (color_row) ), function (x) {unique (color_row [, x]) })
names (gene_list) <- WG$colors2labels (colnames (color_row), prefix='GC')

# We recommend calculating the module scores instead of using eigengene to
# quantify the expression level of particular WGNCA clusters
WG_data <- KG$get_module_score (all_data, paste (sup_save_dir, 
                'WGCNA/Data_module_score.csv', sep='/'), pgenes=gene_list, append_meta=T
)

DEG$seurat_heat (WG_data, cluster_col=T, #get column dendrogram
                 group.by=c('broad_type','date'),
                 normalize=T
)

# ----------cell-cell pairwise correlation----------
# We recommend running correlation on a subtype of data as it is quite
# computationally intensive
PU <- modules::use ('pseudotime_utils.R')
select_cells2 <- all_data$date == 'in_vitro'
select_cells1 <- !select_cells2

# you may choose to run 'partial_corr' or 'distance' to quantify partial
# correlation and Euclidean distance respectively
all_cor <- PU$compute_all_cor (all_data, 
                               method= 'correlation', 
                               assay='RNA',
                               select_cells1=select_cells1
                               select_cells2=select_cells2,
)
# in heatmap
PU$cell_heat (all_cor, all_data@meta.data, 
              features=c('broad_type', 'revised') #group by features
              # first one for select_cells1 (rows), second one for
              # select_cells2 (columns)
)
# in violin plot
PU$cell_violin (all_cor, all_data@meta.data, c('broad_type', 'revised'), normalize=T, num_col=2)

# ----------quadratic programming----------
library (DeconRNASeq)
QP <- modules::use ('quad_prog.R')

select_types <- c('hTSC_OKAE', 'hTSC_TURCO', 'hESC', 'hESC_YAN')
compare_types <- c('ICM', 'TB', 'CTB', 'STB', 'EVT')
cell_sim <- QP$get_cell_similarity (all_data, 
                                    group.by=c('revised', 'broad_type'), 
                                    # in the order of select_types, compare_types
                                    select_types= select_types, 
                                    compare_types=compare_types
)
QP$plot_cell_similarity (all_data, cell_sim, group.by= 'revised', DR='pca')
QP$dim_plot_cell_similarity (all_data, cell_sim, group.by='revised')

# ----------diffusion pseudotime----------
# We did not find diffusion pseudotime suitable for our merged dataset.
# Therefore, we have not written functions for it.
library (destiny)
var_genes <- VariableFeatures (all_data)
GetAssayData (all_data [var_genes, ], assay='RNA', slot='data') %>% as.matrix () -> exp_mat
# The two function below may take a long time to run
dm <- DiffusionMap(t(exp_mat), n_pcs=50)
dpt <- DPT(dm)

dpt@dm@eigenvectors [,1:2] %>% data.frame () %>%
        add_column (all_data$broad_type) %>%
        magrittr::set_colnames (c('DC1', 'DC2', 'Type')) %>%
        ggplot (aes (x=DC1, y=DC2, color=Type))+ geom_point () +
        theme_void () -> plot1
plot2 <-plot(dpt, root = 2, col_by = 'branch')
plot3 <- plot (dpt)
gridExtra::grid.arrange (plot1, plot2, plot3, ncol=2)
