# This script covers dimensionality reduction, DE gene analysis and enrichment
# of GO and KEGG terms
# load all modules and libraries
root_dir <- '/mnt/c/Users/Yutong/Documents/bioinformatics/reproduction/'
setwd (paste (root_dir, 'SingleCellR/utils', sep='/'))
library (tidyverse)
library (Seurat)
library (GOSemSim)

DIR <- modules::use ('dim_red.R')
ML <- modules::use ('marker_list.R')
DEG <- modules::use ('DE_gene.R')
KG <- modules::use ('KEGG_path.R')

# load data
root <- paste (root_dir, 'results/', sep='/')
merge_dir <- paste (root, 'XLYBPZ_Dylan_dir', sep='/')
save_dir <- paste (root, 'manuscript/figure1', sep='/')
sup_save_dir <- paste (root, 'manuscript/figureS1', sep='/')
x <- load (paste (merge_dir, 'final_merged.Robj', sep='/') )
all_data <- get (x)

# plotting the normalised expression of selected genes
plot_genes <- c('POU5F1', 'SOX2', 'NANOG', 'GATA3', 'CDX2', 'TFAP2C')
DEG$seurat_violin (all_data, features=plot_genes, group.by='broad_type')

# ----------dimensionality reductions----------
# This function performs variable gene selection, normalisation, scaling, PCA
# and other dimensionality reduction techniques
all_data <- DIR$run_dim_red (all_data, 
                             normalize=F # whether to use NormalizeData
                             var_scale=T, # whether to scale the dataset
                             run_diff_map=T, # run diffusion map using destiny
                             )

# compare to the Seurat DimPlot function, this function can show multiple
# subplots with different metadata features
highlight <- all_data$date == 'in_vitro'
DIR$plot_dim_red (all_data, 
                  by_group= c('broad_type', 'date', 'paper'), 
                  DR='pca', #can be 'umap' or 'diff_map', default 'pca'
                  dims=c(1,2), #select which dimensionality to show
                  nudge_ratio = 1.2, # optional: adjust the labels of axis relative to the arrow tip
                  num_col = 3, # optional: how many columns to display the subplots
                  size_highlight = highlight # optional: highlight particular cells
)

# This function also works for 3D
DIR$plot_dim_red (all_data, 
                  by_group= c('broad_type', 'date', 'paper'), 
                  DR='pca', #can be 'umap' or 'diff_map'
                  dims=c(1,2,3), #select which dimensionality to show
                  show_axes=T, #optional: whether to show 3D axes
                  all_theta=50, #optional: adjust the azimuthal angle
                  all_phi = 20, #optional: adjust the vertical angle
                  show_label=F, #optional: prevent over-crowding of labels
                  num_col = 3 # optional: how many columns to display the subplots
)

# ----------differentially expressed (DE) genes ----------
# finding DE genes across the entire scRNA-seq dataset can be computationally
# intensive. To use the following function, we make it compulsory to save the
# results.
markers <- DEG$find_DE_genes (all_data, 
                              save_dir, #to save the results
                              feature='broad_type', #DE genes over which feature
                              label='all' # optional: unique tag to the result file
)
# visualise the DE genes in PCA using previous results
# If previous results are not available, this function will call on
# `DEG$find_DE_genes` automatically to generate new results
DIR$plot_gene_PC (all_data, 
                  directory=save_dir, 
                  color_by='broad_type', 
                  label='all', # the label must match the previous DE results
                  show_markers='TF' # optional: show transcriptional factors. 
                  #Turn off this feature with 'show_markers=NULL'
)
# select the top DE genes for visualisation
DE_genes <- DEG$unique_DE_genes (markers, 6)
DE_genes %>% select (group, feature) %>% deframe () -> show_genes
DEG$seurat_heat (all_data, color_row=show_genes, 
                 group.by = c('broad_type', 'date'), #color bars across the rows
                 normalize=T, # row scaling the heatmap
                 col_rotation=90 # rotate the column labels if they are over-crowded
)

# ----------GO/KEGG/Reactome----------
d <- godata('org.Hs.eg.db', ont="BP") #which GO database to use. 
# Here we use 'BP' for biological process, organism is human
# select cells in the trophoblast lineages
markers %>% filter (!group %in% ML$non_TB_lineage) -> tb_markers

# run over-representation test
# Refer to: http://yulab-smu.top/clusterProfiler-book/
# You may also choose 'GO' or 'Reactome' in the 'enrich_area' argument
kk <- KG$compare_cluster_enrichment (tb_markers, d, enrich_area='KEGG')

# visualise in piechart like graph
KG$display_cluster_enrichment (kk, 
                               show_graph='emap', 
                               show_num=20, #optional: how many terms to show
                               default_theme=F #optional: turn off the default ggplot theme
) 

# alternatively, visualise using dotplot
KG$display_cluster_enrichment (kk, show_graph='dotplot', show_num=20) 

# run GSEA for GO/KEGG/Reactome
# again we make it compulsory to save the results as GSEA is quite
# computationally intensive
psea_df <- KG$run_GSEA_all_types (tb_markers, 
                                  enrich_area= 'KEGG',
                                  save_path = paste (sup_save_dir, 
                                  'PSEA_broad_type_all.csv', sep='/')
)

# visualise the results in ridgeplot
KG$ridge_all_types (psea_df, 
                    show_num=40, #optional: how many terms to show
                    not_show_prop=0.01 # optional: remove white space at the ends of the axis
                    # higher values remove more white space
) 
