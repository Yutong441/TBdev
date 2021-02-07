# generate the supplementary figure 1 of the manuscript
# unified atlas of trophoblast development
# all human in vivo datasets, no in vitro cells

setwd('..')
roxygen2::roxygenise()
devtools::load_all('.', export_all=F)
library (tidyverse)
library (org.Hs.eg.db)
library (GOSemSim)

root_dir <- '/mnt/c/Users/Yutong/Documents/bioinformatics/reproduction/'
root <- paste (root_dir, 'results/', sep='/')
merge_dir <- paste (root, 'XLYBPZ_Dylan_dir', sep='/')
save_dir <- paste (root, 'manuscript/figure1', sep='/')
sup_save_dir <- paste (root, 'manuscript/figureS1', sep='/')
all_data <- get (load (paste (merge_dir, 'final_merged_tb.Robj', sep='/') ))

# ----------figure A-C----------
data (CT)
TB_data <- all_data [, !c(all_data$revised %in% CT$in_vitro_cells)]
TB_data$dataset <- gsub ('_[0-9]+$', '', TB_data$paper)
p1 <- plot_dim_red (TB_data, group.by= c('revised', 'date', 'dataset'),
                    DR='pca' , dims=c(1,2), return_sep=T, nudge_ratio=0.3, 
                    plot_type='dim_red_sim', seg_color='black')
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
                                      DR='pca', nudge_ratio=0.3, return_sep=T, plot_type='dim_red_sim')
        paper <- unique (all_datasets[[i]]$paper)
        p2 [[i]] <- one_plot[[1]] + ggtitle (paper) + labs (fill='')
}

p2_final <- gridExtra::grid.arrange (grobs=p2, ncol=3, padding = unit(0.01, "line"))

# ----------figure J----------
all_data2 <- all_data [, ! (all_data$revised %in% c(CT$non_emb_lineage,
                                                    CT$in_vitro_cells))]
data (lineage_markers)
show_genes <- lineage_markers [names (lineage_markers) != 'STR' ]
p3 <- seurat_heat (all_data2, color_row=show_genes, group.by = c('broad_type'), 
                   slot='data', heat_name='norm count',
                   column_legend_labels=c('cell type'),
                   row_legend_labels='lineage markers',
                   column_rotation=90, row_scale=T, center_scale=T,
                   automatic=F)

# ----------figure M and N----------
markers <- find_DE_genes (TB_data, save_dir, group.by='broad_type', label='pairwise', method='pairwise')
tb_ctb <- markers %>% filter (group == 'CTB' & compare_group == 'TB') 
psea_tb_ctb <- run_GSEA_all_types (tb_ctb, org.Hs.eg.db, enrich_area='reactome',
                                   save_path=paste (save_dir, 'RSEA_TB_CTB.csv', sep='/'))
p4 <- enrich_bar (psea_tb_ctb, org.Hs.eg.db, show_num=4, markers=tb_ctb,
            show_gene_labels=6, extend_axis_pos=1.2, extend_axis_neg=4, nudge_x=0.2, shrink_ratio=0.7)

ctb_stb <- markers %>% filter (group == 'STB' & compare_group == 'CTB') 
psea_ctb_stb <- run_GSEA_all_types (ctb_stb, org.Hs.eg.db, enrich_area='reactome',
                                    save_path=paste (save_dir, 'RSEA_CTB_STB.csv', sep='/'))
p5 <- enrich_bar (psea_ctb_stb, org.Hs.eg.db, show_num=4, markers=ctb_stb,
            show_gene_labels=6, extend_axis_pos=1.2, extend_axis_neg=2.2, nudge_x=0.12, shrink_ratio=0.7)


ctb_evt <- markers %>% filter (group == 'EVT' & compare_group == 'CTB') 
psea_ctb_evt <- run_GSEA_all_types (ctb_evt, org.Hs.eg.db, enrich_area='reactome',
                                    save_path=paste (save_dir, 'RSEA_CTB_EVT.csv', sep='/'))
p6 <- enrich_bar (psea_ctb_evt, org.Hs.eg.db, show_num=4, markers=ctb_evt,
            show_gene_labels=6, extend_axis_pos=1.6, extend_axis_neg=2., nudge_x=0.2, shrink_ratio=0.7)

# ----------integration----------
grob_list <- list (p1[[1]]+labs (fill='original \n labels'), p1[[2]], p1[[3]], 
                   p2_final,
                   p3, p4, p5, p6)
lay_mat <- matrix(c(1, 1, 2, 2, 3, 3, 
                    4, 4, 4, 4, 4, 4,
                    4, 4, 4, 4, 4, 4,
                    5, 5, 5, 5, 5, 5,
                    5, 5, 5, 5, 5, 5,
                    6, 6, 7, 7, 8, 8),
                  nrow=6) %>% t()
arrange_plots (grob_list, paste (sup_save_dir, 'final_figureS1.pdf', sep='/'), 
                  lay_mat, plot_width=3, plot_height=6)

save_indiv_plots (grob_list, paste (sup_save_dir, 'figureS1', sep='/'),
                  lay_mat, plot_width=3, plot_height=7
)
