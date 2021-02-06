devtools::load_all('..', export_all=F)
library (tidyverse)
library (GOSemSim)
library (org.Hs.eg.db)

# load data
root_dir <- '/mnt/c/Users/Yutong/Documents/bioinformatics/reproduction/'
root <- paste (root_dir, 'results/', sep='/')
merge_dir <- paste (root, 'XLYBPZ_Dylan_dir', sep='/')
save_dir <- paste (root, 'manuscript/figure2', sep='/')
sup_save_dir <- paste (root, 'manuscript/figureS2', sep='/')
all_data <- get (load (paste (merge_dir, 'final_merged_tb.Robj', sep='/') ))

# ----------figure S2A----------
show_meta <- all_data@meta.data [!(is.na(all_data$MGP_PT) | all_data$broad_type %in% c('PE', 'EPI')), ]
p1 <- pseudo_real_time (show_meta, 'date', 'MGP_PT', 'broad_type', lower_b=0) 

# ----------figure S2C-D: 2-WD ----------
data_dir <- paste (root_dir, 'GPLVM/hier_GP_tf1/', sep='/')
pred_all <- read.csv (paste (data_dir, 'result/prediction_matern_500.csv', sep='/' ))
off_set <- min (all_data$ori_MGP_PT, na.rm=T)
pred_all$x <- pred_all$x - off_set
pred_all$branch <- c('EVT_branch', 'STB_branch')[as.factor (pred_all$branch)]
exp_mat <- read.csv (paste (data_dir, 'data/STREAM_data.csv', sep='/' ), row.names=1)
pseudotime <- read.csv (paste (data_dir, 'result/infer_pt_matern.csv', sep='/' ), row.names=1)
exp_mat %>% as.matrix () %>% t () %>% scale () %>% data.frame () %>% 
        add_column (pseudotime=pseudotime [,'pt_mean'] - off_set) -> exp_mat

markers <- get_DE_from_KL (pred_all, 'EVT_branch', 'STB_branch', divergence='WD')
markers$sign_diver <- sign (markers$logFC)*markers$divergence

show_meta <- all_data@meta.data [!is.na (all_data$MGP_PT) & !all_data$broad_type %in% c('EPI', 'PE') ,]
markers %>% slice_max (sign_diver, n=9) %>% dplyr::select (feature) %>% deframe() -> DE_genes1
p4 <- gene_over_pseudotime (pred_all, exp_mat, DE_genes1, metadata=show_meta, 
                            color_feature = 'broad_type', num_col=3)+ggtitle ('evt_branch')
p4

markers %>% slice_min (sign_diver, n=9) %>% dplyr::select (feature) %>% deframe() -> DE_genes2
p3 <- gene_over_pseudotime (pred_all, exp_mat, DE_genes2, metadata=show_meta, 
                            color_feature = 'broad_type', num_col=3)+ggtitle ('STB_branch')

# ----------figure S3B----------
markers$group [markers$sign_diver <0 ] <- 'STB_branch'
markers$group [markers$sign_diver ==0 ] <- 'none'
markers$group <- partial_relevel (markers$group)
markers$logFC <- markers$divergence
thres <- quantile (markers$logFC, 0.7)

d <- godata(org.Hs.eg.db, ont="BP")
kk <- compare_cluster_enrichment (markers, d, org.Hs.eg.db, enrich_area='KEGG', log_FC_thres=thres)
set.seed(100)
p2 <- display_cluster_enrichment (kk, show_graph='emap', feature_vec=
                                     markers$group, show_num=20) + labs (fill='')

# ----------figure S2E-F----------
#d <- godata(org.Hs.eg.db, ont="BP")
peak_plot <- read.csv (paste (save_dir, 'Switch_STB.csv', sep='/'), row.names=1)
#peak_plot$logFC <- peak_plot$val
#thres <- quantile (peak_plot$logFC, 0.7)
#kk <- compare_cluster_enrichment (peak_plot, d, org.Hs.eg.db, enrich_area='KEGG', log_FC_thres=thres)
#display_cluster_enrichment (kk, show_graph='emap', feature_vec=
#                                     peak_plot$group, show_num=30) + labs (fill='stage')
#MAPK <- gene_per_term (kk, 'MAPK', org.Hs.eg.db, return_val=T)[[1]]
MAPK_sel <- c('GADD45G', 'RAP1B', 'DUSP8', 'RASA1')
#JAK <- gene_per_term (kk, 'JAK-STAT', org.Hs.eg.db, return_val=T)[[1]]
#JAK_sel <- peak_plot [JAK,] %>% slice_max(val, n=4) %>% dplyr::select (feature) %>% deframe ()
JAK_sel <- c('PRLR', 'CSF2RB', 'GH2', 'CSH2')

peak_plot_EVT <- read.csv (paste (save_dir, 'Switch_EVT.csv', sep='/'), row.names=1)
#peak_plot_EVT$logFC <- peak_plot_EVT$val
#thres <- quantile (peak_plot_EVT$logFC, 0.7)
#kk_EVT <- compare_cluster_enrichment (peak_plot_EVT, d, org.Hs.eg.db, enrich_area='KEGG', log_FC_thres=thres)
#display_cluster_enrichment (kk_EVT, show_graph='emap', feature_vec=
#                                     peak_plot_EVT$group, show_num=30) + labs (fill='stage')
#
#PI3K <- gene_per_term (kk_EVT, 'PI3K-Akt', org.Hs.eg.db, return_val=T)[[1]]
#PI3K_sel <- peak_plot_EVT [PI3K,] %>% slice_max(val, n=4) %>% dplyr::select (feature) %>% deframe ()
PI3K_sel <- c('FLT4', 'LAMA4', 'CSF1R', 'EFNA1')

# obtain expresion matrix
p5 <- gene_over_pseudotime (pred_all, exp_mat, MAPK_sel, show_meta,
                            color_feature = 'broad_type', num_col=2,
                            peak_data=peak_plot)+ggtitle ('MAPK')
p6 <- gene_over_pseudotime (pred_all, exp_mat, JAK_sel, show_meta,
                            color_feature = 'broad_type', num_col=2,
                            peak_data=peak_plot)+ggtitle ('JAK-STAT')
p7 <- gene_over_pseudotime (pred_all, exp_mat, PI3K_sel, show_meta,
                            color_feature = 'broad_type', num_col=2,
                            peak_data=peak_plot_EVT)+ggtitle ('PI3K-AKT')

# ----------figure S2E----------
save_dir4 <- paste (root, 'manuscript/figure4', sep='/')
color_row <- read.csv ( paste ( save_dir4, 
                        'WGCNA/module_genes.csv' , sep='/'), row.names=1)
gene_list <- lapply (as.list (colnames (color_row) ), function (x) {
                             unique (color_row [, x]) })
ori_names <- colnames (color_row)
names (gene_list) <- colors2labels (colnames (color_row), prefix='GC')

invivo <- all_data [, all_data$date != 'in_vitro']
module_score <- get_module_score (all_data, append_meta=T, paste (save_dir4, 
                'WGCNA/Data_module_score.csv', sep='/'), pgenes=gene_list)

vivo_mod <- module_score [, module_score$date != 'in_vitro']
vivo_scale <- scale_seurat (vivo_mod, row_scale=T, slot_data='counts')
module_data <- vivo_scale [, !is.na (vivo_scale$MGP_PT) & !vivo_scale$broad_type %in% c('EPI', 'PE')]
sel_genes <- rownames (module_data)
names (sel_genes) <- c('pre-implant', 'non-specific', 'STB', 'non-specific', 'EPI', 'CTB', 'ICM', 'ICM',
                       'non-specific', 'EVT', 'cleavage')
row_levels <- partial_relevel (names (sel_genes)) %>% levels()
p8 <- seurat_heat (module_data, group.by=c('epil_branch','broad_type'),
                 row_scale=T, color_row= sel_genes,
                 column_reorder_levels = list (format_conf$branch_order, format_conf$cell_order),
                 column_legend_labels= c('branch', 'cell type'), 
                 show_column_bars=c(F,T),
                 row_legend_labels='WGCNA clusters',
                 cluster_rows=T, heat_name='norm count', center_scale=T,
                 group_order = order (module_data$MGP_PT),
                 automatic=F, left_HA=F, slot_data='counts'
)

grob_list <- list (p1, p2, p3, p4, p5, p6, p7, p8)
lay_mat <- matrix(c(1, 1, 1, 2, 2, 2, 
                    3, 3, 3, 4, 4, 4,
                    5, 5, 6, 6, 7, 7,
                    8, 8, 8, 8, 8, 8
                    ),
                  nrow=6) %>% t()
arrange_plots (grob_list, paste (sup_save_dir, 'final_figureS2.pdf', sep='/'), 
                  lay_mat, plot_width=3, plot_height=7)

# find peak genes
save_peak_genes (peak_plot, paste (sup_save_dir, 'PT_marker/Peak_EVT_PI3K.csv', sep='/'), 'PI3K-Akt', kk, max_num=50)
save_peak_genes (peak_plot, paste (sup_save_dir, 'PT_marker/Peak_EVT_HIF1.csv', sep='/'), 'HIF-1', kk, max_num=50)
save_peak_genes (peak_plot, paste (sup_save_dir, 'PT_marker/Peak_EVT_relaxin.csv', sep='/'), 'Relaxin', kk, max_num=50)
save_peak_genes (peak_plot, paste (sup_save_dir, 'PT_marker/Peak_EVT_all.csv', sep='/'))
