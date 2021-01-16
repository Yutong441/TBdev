# generate the supplementary figure 1 of the manuscript
# unified atlas of trophoblast development
# all human in vivo datasets, no in vitro cells

setwd('..')
roxygen2::roxygenise()
devtools::load_all()
library (tidyverse)
library (org.Hs.eg.db)
library (GOSemSim)

root_dir <- '/mnt/c/Users/Yutong/Documents/bioinformatics/reproduction/'
root <- paste (root_dir, 'results/', sep='/')
merge_dir <- paste (root, 'XLYBPZ_Dylan_dir', sep='/')
save_dir <- paste (root, 'manuscript/figure1', sep='/')
sup_save_dir <- paste (root, 'manuscript/figureS1', sep='/')
x <- load (paste (merge_dir, 'final_merged_tb.Robj', sep='/') )
all_data <- get (x)

# ----------figure A-C----------
data (CT)
TB_data <- all_data [, !c(all_data$revised %in% CT$in_vitro_cells)]
TB_data$dataset <- gsub ('_[0-9]+$', '', TB_data$paper)
p1 <- plot_dim_red (TB_data, group.by= c('revised', 'date', 'dataset'),
                    DR='pca' , dims=c(1,2), return_sep=T, nudge_ratio=0.3)
# ----------figure D-I----------
save_robj <- c('Xiang_2019/Xiang_R.Robj', 'Liu_2018/Liu_R.Robj',
               'Blakeley_2015/Blakeley_R.Robj',
               'Petropoulos_2016/Petropoulos_R.Robj', 'Yan_2013/Yan_R.Robj',
               'Zhou_2019/Zhou_R.Robj')

merge_dir <- sapply (save_robj, function(x){strsplit(x, '')[[1]][[1]]})
data_dir <- paste (root_dir, 'data/', sep='/')
merge_dir <- paste (c(data_dir, merge_dir, '_dir'), collapse='')
if (!dir.exists (merge_dir)) {dir.create (merge_dir)}
all_datasets <- load_all_data (save_robj, data_dir)

# run PCA on each dataset
for (i in 1:length (all_datasets)){
        all_datasets [[i]] <- clean_metadata (all_datasets[[i]])
        all_datasets[[i]] <- run_dim_red (all_datasets[[i]], run_diff_map=F, var_scale=T,
                             normalize=F, find_var_features=T, run_umap=F)
}
p2 <- list ()
for (i in 1:length (all_datasets)){
        one_plot <- plot_dim_red (all_datasets[[i]], group.by= c('Type'),
                                      DR='pca', all_labels=T, nudge_ratio=0.3, return_sep=T)
        paper <- unique (all_datasets[[i]]$paper)
        p2 [[i]] <- one_plot[[1]] + ggtitle (paper) + labs (fill='')
}

gridExtra::grid.arrange (grobs=p2, ncol=3)

# ----------figure J----------
all_data2 <- all_data [, ! (all_data$revised %in% c(CT$non_emb_lineage,
                                                    CT$in_vitro_cells))]
data (lineage_markers)
show_genes <- lineage_markers [names (lineage_markers) != 'STR' ]
p3 <- seurat_heat (all_data2, color_row=show_genes, group.by = c('revised'), 
                   slot='data', heat_name='norm count',
                   column_legend_labels=c('cell type'),
                   column_rotation=90, row_scale=T, center_scale=T)
p3

# ----------figure M and N----------
markers <- find_DE_genes (all_data, save_dir, feature='broad_type', label='all')
markers %>% filter (! group %in% CT$non_TB_lineage ) -> tb_markers

d <- godata(org.Hs.eg.db, ont="BP")
ra <- compare_cluster_enrichment (tb_markers, d, org.Hs.eg.db, enrich_area='reactome')
p6 <- display_cluster_enrichment (ra, show_graph='emap', feature_vec=tb_markers$group, 
                                     show_num=20, vert_just=2.5) + labs (fill = '')
p6

rsea_df <- run_GSEA_all_types (tb_markers, enrich_area='reactome', 
                                  save_path = paste (sup_save_dir, 'RSEA_broad_type_vivo.csv', sep='/'))
p7<- ridge_all_types (rsea_df, show_num=20)+ labs (fill='')
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
arrange_plots (grob_list, paste (sup_save_dir, 'final_figureS1.pdf', sep='/'), 
                  lay_mat, plot_width=3, plot_height=6)
