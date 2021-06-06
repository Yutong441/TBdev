# generate the supplementary figure 1 of the manuscript
# unified atlas of trophoblast development
# all human in vivo datasets, no in vitro cells

setwd('..')
#roxygen2::roxygenise()
devtools::load_all('.', export_all=F)
library (tidyverse)
library (org.Hs.eg.db)

root_dir <- '/mnt/c/Users/Yutong/Documents/bioinformatics/reproduction/'
root <- paste (root_dir, 'results/', sep='/')
merge_dir <- paste (root, 'XLYBPZ_Dylan_dir', sep='/')
save_dir <- paste (root, 'manuscript/figure1', sep='/')
sup_save_dir <- paste (root, 'manuscript/figureS1', sep='/')
all_data <- get (load (paste (merge_dir, 'final_merged_tb.Robj', sep='/') ))
TB_data <- all_data [, !c(all_data$revised %in% CT$in_vitro_cells)]

# ----------figure A-C----------
data (CT)
TB_data$dataset <- gsub ('_[0-9]+$', '', TB_data$paper)
p1 <- plot_dim_red (TB_data, group.by= c('revised', 'date', 'dataset'),
                    DR='pca' , dims=c(1,2), return_sep=T, nudge_ratio=0.3, 
                    plot_type='dim_red_sim', seg_color='black',
                    nudge_dimname=0.2, nudge_ortho=0.3)

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
                                  DR='pca', nudge_ratio=0.3, return_sep=T, plot_type='dim_red_sim',
                                  nudge_dimname=0.2, nudge_ortho=0.3)
        paper <- unique (all_datasets[[i]]$paper)
        p2 [[i]] <- one_plot[[1]] + ggtitle (paper) + labs (fill='')
}

p2_final <- gridExtra::grid.arrange (grobs=p2, ncol=3, padding = unit(0.01, "line"))

# ----------figure J----------
all_data2 <- all_data [, ! (all_data$revised %in% c(CT$non_emb_lineage,
                                                    CT$in_vitro_cells))]
data (lineage_markers)
show_genes <- lineage_markers [names (lineage_markers) != 'STR' ]
show_genes <- c(show_genes, AM='ISL1', AM='GABRP', AM='VTCN1', 
                HLA='HLA-A', HLA='HLA-B', HLA='HLA-C')

data (format_conf)
new_color <- list (color_vec=c(format_conf$color_vec, 'AM'='#0352fc', 'HLA'='#5f6675'))

p3 <- seurat_heat (all_data2, color_row=show_genes, group.by = c('broad_type'), 
                   slot='data', heat_name='norm count',
                   column_legend_labels=c('cell type'),
                   row_legend_labels='lineage markers',
                   column_rotation=90, row_scale=T, center_scale=T,
                   automatic=F, AP=new_color)

# ----------figure M and N----------
# volcano plot of TF
data (TF)
markers <- find_DE_genes (TB_data, save_dir, group.by='broad_type',
                          label='pairwise', method='pairwise')
tf_mark <- markers %>% filter (feature %in% TF)

p4<- seurat_volcano (tf_mark, 'CTB', 'TB', weighting='logFC',
                nudge_y=80, length_ratio=0.8)
p5<- seurat_volcano (tf_mark, 'STB', 'CTB', weighting='logFC',
                nudge_y=80, length_ratio=0.8, nudge_x=0.1)
p6<- seurat_volcano (tf_mark, 'EVT', 'CTB', weighting='logFC',
                nudge_y=80, length_ratio=0.8, nudge_x=0.1)

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
