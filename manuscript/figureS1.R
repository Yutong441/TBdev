# generate the supplementary figure 1 of the manuscript
# unified atlas of trophoblast development
# all human in vivo datasets, no in vitro cells

root_dir <- '/mnt/c/Users/Yutong/Documents/bioinformatics/reproduction/'
setwd (paste (root_dir, 'SingleCellR/utils', sep='/'))
library (modules)
library (tidyverse)
library (Seurat)
library (GOSemSim)
DIR <- modules::use ('dim_red.R')
BU <- modules::use ('batch_utils.R')
ML <- modules::use ('marker_list.R')
DEG <- modules::use ('DE_gene.R')
PD <- modules::use ('plot_DR.R')
CC <- modules::use ('cycle_analysis.R')
KG <- modules::use ('KEGG_path.R')
FO <- modules::use ('format.R')

root <- paste (root_dir, 'results/', sep='/')
merge_dir <- paste (root, 'XLYBPZ_Dylan_dir', sep='/')
sup_save_dir <- paste (root, 'manuscript/figureS1', sep='/')
x <- load (paste (merge_dir, 'final_merged_vivo.Robj', sep='/') )
all_data <- get (x)

# select in vivo cells, keep all other cells
#all_data <- all_data [, ! (all_data$revised %in% c(ML$non_emb_lineage))]
#all_data <- DIR$run_dim_red (all_data, run_diff_map=F, var_scale=T,
#                             normalize=F, find_var_features=T, run_umap=F)
#all_data <- ML$clean_metadata (all_data)

# ----------figure A-C----------
TB_data <- all_data [, !c(all_data$revised %in% ML$in_vitro_cells)]
TB_data$dataset <- gsub ('_[0-9]+$', '', TB_data$paper)
p1 <- PD$gg_plot_dim_red (TB_data, by_group = c('revised', 'date', 'dataset'),
                          DR='pca' , dims=c(1,2), return_sep=T, nudge_ratio=1.5)

# ----------figure D-I----------
save_robj <- c('Xiang_2019/Xiang_R.Robj', 'Liu_2018/Liu_R.Robj',
               'Blakeley_2015/Blakeley_R.Robj',
               'Petropoulos_2016/Petropoulos_R.Robj', 'Yan_2013/Yan_R.Robj',
               'Zhou_2019/Zhou_R.Robj')

merge_dir <- sapply (save_robj, function(x){strsplit(x, '')[[1]][[1]]})
data_dir <- paste (root_dir, 'data/', sep='/')
merge_dir <- paste (c(data_dir, merge_dir, '_dir'), collapse='')
if (!dir.exists (merge_dir)) {dir.create (merge_dir)}
all_datasets <- BU$load_all_data (save_robj, data_dir)

# run PCA on each dataset
for (i in 1:length (all_datasets)){
        all_datasets [[i]] <- ML$clean_metadata (all_datasets[[i]])
        all_datasets[[i]] <- DIR$run_dim_red (all_datasets[[i]], run_diff_map=F, var_scale=T,
                             normalize=F, find_var_features=T, run_umap=F)
}
p2 <- list ()
for (i in 1:length (all_datasets)){
        one_plot <- DIR$plot_dim_red (all_datasets[[i]], by_group= c('Type'),
                                      DR='pca', all_labels=T, nudge_ratio=1.5, return_sep=T)
        paper <- unique (all_datasets[[i]]$paper)
        p2 [[i]] <- one_plot[[1]] + ggtitle (paper) + labs (fill='')
}

gridExtra::grid.arrange (grobs=p2, ncol=3)

# ----------figure J----------
all_data2 <- all_data [, ! (all_data$revised %in% c(ML$non_emb_lineage,
                                                    ML$in_vitro_cells))]
show_genes <- ML$lineage_markers [names (ML$lineage_markers) != 'STR' ]
p3 <- DEG$seurat_heat (all_data2, color_row=show_genes, group.by = c('revised'), 
                       slot='data', heat_name='norm count',
                       col_legend_labels=c('cell type'),
                       col_rotation=90, normalize=T, center_scale=T)
p3

# ----------figure M and N----------
save_dir <- paste (root, 'manuscript/figure1', sep='/')
markers <- DEG$find_DE_genes (all_data, save_dir, feature='broad_type', label='all')
markers %>% filter (! group %in% ML$non_TB_lineage ) -> tb_markers

d <- godata('org.Hs.eg.db', ont="BP")
ra <- KG$compare_cluster_enrichment (tb_markers, d, enrich_area='reactome')
p6 <- KG$display_cluster_enrichment (ra, show_graph='emap', feature_vec=tb_markers$group, 
                                     show_num=20, vert_just=2.5) + labs (fill = '')
p6

rsea_df <- KG$run_GSEA_all_types (tb_markers, enrich_area='reactome', 
                                  save_path = paste (sup_save_dir, 'RSEA_broad_type_vivo.csv', sep='/'))
p7<- KG$ridge_all_types (rsea_df, show_num=20)+ labs (fill='')
p7

# ----------integration----------
grob_list <- list (p1[[1]]+labs (fill=''), p1[[2]], p1[[3]], 
                   p2[[1]], p2[[2]], p2[[3]], p2[[4]], p2[[5]], p2[[6]],
                   p3, p6, p7+theme (aspect.ratio=1.5))
lay_mat <- matrix(c(1, 1, 2, 2, 3, 3, 
                    4, 4, 5, 5, 6, 6,
                    7, 7, 8, 8, 9, 9,
                    10,10,10,10,10,10,
                    10,10,10,10,10,10,
                    11,11,11,12,12,12),
                  nrow=6) %>% t()
FO$arrange_plots (grob_list, paste (sup_save_dir, 'final_figureS1.pdf', sep='/'), 
                  lay_mat, plot_width=3, plot_height=6)
