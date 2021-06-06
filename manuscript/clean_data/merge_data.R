# This script merges data for downstream analysis
setwd ('utils')
library (modules)
library (Seurat)
library (tidyverse)
library (topGO)
library (org.Hs.eg.db)
DIR <- modules::use ('dim_red.R')
BU <- modules::use ('batch_utils.R')
ML <- modules::use ('marker_list.R')
EM <- modules::use ('expr_mat.R')
GO <- modules::use ('GO_analysis.R')
DEG <- modules::use ('DE_gene.R')

root <- 'data/'
save_dir <- 'result/'

save_robj <- c('Xiang_2019/Xiang_R.Robj', 'Liu_2018/Liu_R.Robj',
               'Blakeley_2015/Blakeley_R.Robj',
               'Petropoulos_2016/Petropoulos_R.Robj', 'Yan_2013/Yan_R.Robj',
               'Zhou_2019/Zhou_R.Robj')

merge_dir <- sapply (save_robj, function(x){strsplit(x, '')[[1]][[1]]})
merge_dir <- paste (c(save_dir, merge_dir, '_dir'), collapse='')
if (!dir.exists (merge_dir)) {dir.create (merge_dir)}
all_datasets <- BU$load_all_data (save_robj, root)
load (paste (merge_dir, 'merged.Robj', sep='/'))

# ----------unmerged data----------#
unmerged <- BU$merge_seurat (all_datasets, c('RNA'))
unmerged <- DIR$run_dim_red (unmerged, run_diff_map=T, var_scale=T, normalize=F)
unmerged <- ML$clean_metadata (unmerged)

DIR$plot_dim_red (unmerged, by_group = c('paper', 'revised_date'), DR='pca', all_labels=F )
DIR$plot_dim_red (unmerged, by_group = c('paper', 'revised_date'), DR='pca', all_labels=F, dims=c(1,3))
DIR$plot_dim_red (unmerged, by_group= c('revised', 'date'), DR='pca')
DIR$plot_dim_red (unmerged, by_group= c('revised', 'date'), DR='pca', dims=c(2,3))
DIR$plot_dim_red (unmerged, by_group = c('paper', 'revised'), DR='diff_map')
DIR$plot_dim_red (unmerged, by_group= c('revised', 'date'), DR='diff_map')
save (unmerged, file=paste (merge_dir, 'merged.Robj', sep='/'))

# heatmap
DEG$seurat_heat (unmerged, color_row=ML$lineage_markers, group.by = c('revised', 'date' ), slot='data')
DEG$seurat_heat (unmerged, color_row=ML$lineage_markers, group.by = c('date', 'revised' ), slot='data')
DEG$seurat_heat (unmerged, color_row=ML$FGF, group.by = c('revised', 'date' ), slot='data')
DEG$seurat_heat (unmerged, color_row=ML$FGF, group.by = c('date', 'revised' ), slot='data')

DEG$seurat_heat (unmerged, color_row=ML$lineage_markers, group.by = c('seurat_clusters', 'revised'), slot='data')
DEG$seurat_heat (unmerged, color_row=ML$lineage_markers, group.by = c('revised', 'seurat_clusters'), slot='data')
DEG$seurat_heat (unmerged, color_row=ML$lineage_markers, group.by = c('date', 'revised'), slot='data')

DE_genes <- DEG$DE_lineage_genes (unmerged, merge_dir, top_number=3)
DEG$seurat_heat (TB_data, color_row=DE_genes, group.by = c('date', 'revised' ), slot='data')
DEG$seurat_heat (TB_data, color_row=DE_genes, group.by = c('revised', 'date' ), slot='data')

markers <- DEG$find_DE_genes (unmerged, merge_dir, feature='reassign', label='all')
GO$GO_output (markers, rownames (unmerged), merge_dir)

# ----------remove stromal cells----------#
unmerged <- DIR$run_dim_red (unmerged, select_cells = unmerged$Type != 'STR',
                             run_diff_map=F, var_scale=T, run_umap=F)
unmerged <- ML$clean_metadata (unmerged)
DIR$plot_dim_red (unmerged, by_group= c('paper', 'revised_date'), DR='pca', all_labels=F)
DIR$plot_dim_red (unmerged, by_group= c('revised', 'date'), DR='pca')
DIR$plot_dim_red (unmerged, by_group= c('revised', 'date'), DR='pca', dim=c(1,3))
DIR$plot_gene_PC (unmerged, directory= merge_dir, color_by='revised')

non_STR <- DIR$RunDiffusion  (unmerged [, unmerged$Type != 'STR'])
DIR$DimPlot_3D (non_STR, 'Type', DR='diff_map')
DIR$plot_dim_red (non_STR, by_group= c('revised', 'date'), DR='diff_map', dims=c(1, 3))
DIR$plot_gene_PC (non_STR, dims=c(1,3), directory=merge_dir)

# ---------- TB lineage only ----------#
TB_data <- unmerged [, !(unmerged$revised %in% ML$non_TB_lineage) ]
TB_data <- TB_data [, (as.numeric (TB_data$date) >=5 )]
TB_data <- DIR$run_dim_red (TB_data, run_diff_map=F, find_var_features=T)
TB_data <- ML$clean_metadata (TB_data)

DIR$plot_dim_red (TB_data, by_group= c('paper', 'revised_date'), DR='pca', all_labels=F, dims=c(1,2))
DIR$plot_dim_red (TB_data, by_group= c('revised', 'date'), DR='pca', dims=c(1,2))
DIR$plot_dim_red (TB_data, by_group= c('revised', 'date'), DR='pca', dims=c(1,3))
DIR$plot_dim_red (TB_data, by_group= c('revised', 'seurat_clusters'), DR='pca', dims=c(1:2))
DIR$plot_dim_red (TB_data, by_group= c('paper', 'revised'), DR='diff_map')
DIR$plot_dim_red (TB_data, by_group= c('revised', 'date'), DR='diff_map', dims=c(2:3))
DIR$plot_gene_PC (TB_data, dims=c(1,2), directory=merge_dir, color_by='revised')
DIR$plot_gene_PC (TB_data, dims=c(1,2), directory=merge_dir, show_markers='TF', color_by='revised')
table (TB_data$revised, TB_data$seurat_clusters) %>% write.csv (paste (merge_dir, 'celltype_table.csv', sep='/') )

DE_genes <- DEG$DE_lineage_genes (TB_data, merge_dir, top_number=3, label='TB')
DEG$seurat_heat (TB_data, color_row=DE_genes, group.by = c('revised', 'date' ), slot='data')
DEG$seurat_heat (TB_data, color_row=DE_genes, group.by = c('date', 'revised' ), slot='data')
DEG$seurat_heat (TB_data, color_row=DE_genes, group.by = c('seurat_clusters', 'revised' ), slot='data')

# save to csv for GPLVM in python
select_cells <- TB_data$date == 'DNA' | TB_data$revised == 'STR'
EM$save_to_csv (TB_data, merge_dir, assay='RNA', select_cells = !select_cells,
                select_genes=VariableFeatures(unmerged))

# ----------Combine blastoid data with the in vivo data----------
devtools::load_all('../..', export_all=F)
root_dir <- '/mnt/c/Users/Yutong/Documents/bioinformatics/reproduction/'
root <- paste (root_dir, 'results/', sep='/')
merge_dir <- paste (root, 'XLYBPZ_Dylan_dir', sep='/')

blastoid <- readRDS(paste (root_dir, 'results/Yu_2021/blastoid.rds', sep='/') )
blastoid$date <- 'in_vitro'
blastoid$broad_type <- blastoid$Type
blastoid$final_cluster <- blastoid$Type
all_data_a <- get(load(file=paste (merge_dir, 'final_merged_tb.Robj', sep='/')))

Yana <- readRDS(paste (root_dir, 'results/Yanagida_2021/Yanagida_R.rds', sep='/') )
Yana <- Yana [, !Yana$broad_type %in% c('bTB', 'bEPI', 'bPE', 'bEPI-PE')]
all_data <- TBdev::merge_seurat (list (all_data_a, blastoid), assays='RNA')
all_data <- TBdev::merge_seurat (list (all_data, Yana), assays='RNA')
all_data <- all_data [, !all_data$broad_type %in% c('bESC', 'hTSC-TURCO')]

table (all_data$broad_type)
in_vitro_cells <- all_data$date == 'in_vitro'
VariableFeatures (all_data) <- VariableFeatures (all_data_a)
all_data <- run_dim_red (all_data, run_diff_map=F, var_scale=T,cluster=F,
                         normalize=F, find_var_features=F, run_umap=F,
                         select_cells=!in_vitro_cells) #PCA projection for in vitro
all_data <- all_data [, all_data$final_cluster != 'uCTB']

rm (all_data_a)
saveRDS (all_data, file=paste (merge_dir, 'merged_blastoid.rds', sep='/'))
