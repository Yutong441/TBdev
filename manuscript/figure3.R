devtools::load_all('..', export_all=F)
library (tidyverse)

root_dir <- '/mnt/c/Users/Yutong/Documents/bioinformatics/reproduction/'
root <- paste (root_dir, 'results/', sep='/')
merge_dir <- paste (root, 'XLYBPZ_Dylan_dir', sep='/')
save_dir <- paste (root, 'manuscript/figure3', sep='/')
sup_save_dir <- paste (root, 'manuscript/figureS3', sep='/')
all_data_a <- get (load (paste (merge_dir, 'final_merged_tb.Robj', sep='/') ))
all_data <- all_data_a [, all_data_a$assigned_cluster != 'uCTB']

# ----------figure 3A ----------
in_vitro <- all_data [, all_data$date == 'in_vitro' ]
data (lineage_markers)
plot_genes <- lineage_markers [names (lineage_markers)  == 'TB']
plot_genes <- c('POU5F1', 'SOX2', 'NANOG', 'KRT7', 'KRT18', 'TFAP2C',
                'TFAP2A','GATA3','GATA2', 'ELF5', 'CDX2', 'TEAD4', 'HLA-A',
                'HLA-B','TP63')

p1 <- seurat_violin (in_vitro, features=plot_genes, group.by='final_cluster',
                     num_col=5, box_plot=F, free_xy='free_y', lower_b=0)

# ----------figure 3B ----------
TSC_cells <- all_data$revised %in% c('hTSC_OKAE', 'hTSC_TURCO') 
data (CT)
highlight <- all_data$assigned_cluster %in% CT$in_vitro_cells
p2 <- plot_dim_red (all_data, group.by= c('broad_type'), DR='pca' , return_sep=T,
                    size_highlight=highlight, highlight_ratio=2, dims=c(1,2),
                    nudge_ratio=0.2, plot_type='dim_red_sim', seg_color='black', 
                    AP=list (point_edge_color='gray'))

# ----------figure 3C----------
sup_save_dir2 <- paste (root, 'manuscript/figureS2', sep='/')
epg <- read.csv (paste (sup_save_dir2, 'result/STREAM_graph.csv', sep='/'))
show_meta <- all_data@meta.data
show_meta %>% filter (!is.na (MGP_PT) ) %>% filter (!broad_type %in% c('PE', 'EPI') ) -> show_meta

new_name <- c('TB_stem', 'EVT_branch', 'STB_branch', NA)
old_name <- c('S2,S1', 'S0,S1', 'S3,S1', 'blank')
data (format_conf)
epg$branch_name <- partial_relevel (new_name[match (epg$branch, old_name )], format_conf$branch_order)

epg %>% arrange (branch_name) %>% mutate( index = 1:nrow (epg)) %>% 
        group_by (branch) %>% slice_min (index, n=1) %>% drop_na() %>% 
        as.data.frame () %>% dplyr::select (x, y, z, branch_name) %>% 
        magrittr::set_colnames (c('x', 'y', 'z', 'feature')) -> label_epg

show_meta %>% slice_min (PT3, n=nrow(show_meta)-3) %>% 
        normalize_prob(c('hTSC_OKAE', 'hTSC_TURCO')) %>%
        dim_red_3D_traj ('PT1', 'PT2', 'PT3', c('hTSC_OKAE', 'hTSC_TURCO'), epg, 'x', 'y',
        'z', 'branch_name', all_theta=50, all_phi=0, show_label=T, further_repel=T,
        repel_force=0.2, lab_just=c(0.08, 0.02, 0.02), label_col='broad_type',
        label_traj_text=label_epg, num_col=1, AP=list (point_edge_color='white'), 
        hor_just=0.1, dim_vjust=3) + 
        labs (fill='relative probability') +
        scale_color_manual (values=rep('black',3))+guides(color=F) -> p3

# ----------figure 3D: cell-cell correlation ----------
select_cells2 <- all_data$revised %in% CT$in_vitro_cells
select_cells1 <- !(all_data$revised %in% CT$in_vitro_cells)
all_cor <- compute_all_cor (all_data, method= 'correlation', assay='RNA',
                            select_cells2=select_cells2,
                            select_cells1=select_cells1)
p4 <- cell_violin (all_cor, all_data@meta.data, c('final_cluster', 'final_cluster'), 
                   box_plot=T, num_col=2, legend_col=1)

# ----------figure 3E: TF activities----------
save_dir4 <- paste (root, 'manuscript/figure4', sep='/')
color_row <- read.csv ( paste ( save_dir4, 
                        'WGCNA/module_genes.csv' , sep='/'), row.names=1)
gene_list <- lapply (as.list (colnames (color_row) ), function (x) {
                             unique (color_row [, x]) })
ori_names <- colnames (color_row)
names (gene_list) <- colors2labels (colnames (color_row), prefix='GC')

module_score <- get_module_score (all_data, append_meta=T, paste (save_dir4, 
                'WGCNA/Data_module_score.csv', sep='/'), pgenes=gene_list)
vitro_mod <- module_score [, module_score$date=='in_vitro']
sel_genes <- rownames (module_score)
names (sel_genes) <- c('pre-implant', 'non-specific', 'STB', 'non-specific', 'EPI', 'CTB', 'ICM', 'ICM',
                       'non-specific', 'EVT', 'cleavage')
row_levels <- partial_relevel (names (sel_genes)) %>% levels()
p5 <- seurat_heat (vitro_mod, group.by=c('broad_type'),
                 row_scale=T, color_row= sel_genes,
                 row_legend_labels='WGCNA clusters',
                 column_legend_labels='cell type',
                 cluster_rows=T, heat_name='norm count', center_scale=T,
                 automatic=F, left_HA=F, slot_data='counts',
                 column_title_fontface='plain',
                 column_title_side='bottom',
                 main_width=6, main_height=14, column_rotation=90
)

# ----------figure 3F: pathway module----------
sup_save_dir2 <- paste (root, 'manuscript/figureS2', sep='/')
module_scores <- get_module_score (all_data, save_path=paste (sup_save_dir2, 
                                        'Data_module_scores.csv', sep='/'))
colnames (module_scores) <- gsub ('^X', '', colnames (module_scores))
rownames (module_scores) <- gsub ('\\.', '-', rownames (module_scores) )
meta_data <- all_data_a@meta.data [match (colnames (module_scores), colnames (all_data_a) ), ]
module_signal <- Seurat::CreateSeuratObject ( module_scores, meta.data = meta_data )
module_signal <- module_signal [, !module_signal$broad_type %in% c('EPI', 'PE', 'hESC', 'hESC-YAN')]
module_signal <- module_signal [, !module_signal$final_cluster %in% c('uCTB')]

sel_seurat <- average_by_group (module_signal, 'final_cluster', rownames (module_signal))
p6<-seurat_heat (sel_seurat, 'final_cluster', rownames (module_signal),
             main_width=7, main_height=13, column_split=NA,
             column_rotation=90, show_column_names=T, cluster_column=T,
             center_scale=T, column_legend_labels=c('cell type'), row_scale=T,
             grid_height=5, heat_grid_height=6, automatic=F
)

# arrange all figures
grob_list <- list (p1 + theme (aspect.ratio=1., 
                        axis.title.x=element_blank(), legend.position='top') + labs (fill =''), 
                   p2[[1]] +labs (fill=''), 
                   p3+ theme(legend.position='top'), 
                   p4+ labs (fill='') + theme (axis.title.x=element_blank(), aspect.ratio=0.5), 
                   p5, p6)
lay_mat <- matrix(c(1, 1, 1, 2, 2, 2, 
                    3, 3, 4, 4, 4, 4,
                    3 ,3 ,5 ,5 ,6 ,6
                    ),
                  nrow=6) %>% t()
arrange_plots (grob_list, paste (save_dir, 'final_figure3.pdf', sep='/'),
               lay_mat, plot_width=3.5, plot_height=8)

save_indiv_plots (grob_list, paste (save_dir, 'figure3', sep='/'),
               lay_mat, plot_width=3.5, plot_height=8)
