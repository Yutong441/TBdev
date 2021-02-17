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
all_data$realtime <- gsub ('^D', 'E', all_data$date)
all_data$realtime <- partial_relevel (all_data$realtime)
show_meta <- all_data@meta.data [!(is.na(all_data$MGP_PT) | all_data$broad_type %in% c('PE', 'EPI')), ]
p1 <- pseudo_real_time (show_meta, 'realtime', 'MGP_PT', 'broad_type', lower_b=0) 

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
markers %>% slice_max (sign_diver, n=10) %>% dplyr::select (feature) %>% deframe() -> DE_genes1
p4 <- gene_over_pseudotime (pred_all, exp_mat, DE_genes1, metadata=show_meta, 
                            color_feature = 'broad_type', num_col=2)+ggtitle ('EVT_branch')

markers %>% slice_min (sign_diver, n=10) %>% dplyr::select (feature) %>% deframe() -> DE_genes2
p3 <- gene_over_pseudotime (pred_all, exp_mat, DE_genes2, metadata=show_meta, 
                            color_feature = 'broad_type', num_col=2)+ggtitle ('STB_branch')

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
d <- godata(org.Hs.eg.db, ont="BP")
peak_plot <- read.csv (paste (save_dir, 'Switch_STB.csv', sep='/'), row.names=1)
peak_plot$logFC <- peak_plot$val
thres <- quantile (peak_plot$logFC, 0.7)
#kk <- compare_cluster_enrichment (peak_plot, d, org.Hs.eg.db, enrich_area='KEGG', log_FC_thres=thres)
#display_cluster_enrichment (kk, show_graph='emap', feature_vec=
#                                     peak_plot_EVT$group, show_num=30) + labs (fill='stage')

#MAPK <- gene_per_term (kk, 'MAPK', org.Hs.eg.db, return_val=T)[[1]]
#JAK <- gene_per_term (kk, 'JAK-STAT', org.Hs.eg.db, return_val=T)[[1]]
#JAK_sel <- peak_plot [JAK,] %>% slice_max(val, n=4) %>% dplyr::select (feature) %>% deframe ()
JAK_sel <- c('PRLR', 'CSF2RB', 'GH2', 'CSH2')

peak_plot_EVT <- read.csv (paste (save_dir, 'Switch_EVT.csv', sep='/'), row.names=1)
peak_plot_EVT$logFC <- peak_plot_EVT$val
thres <- quantile (peak_plot_EVT$logFC, 0.7)
kk_EVT <- compare_cluster_enrichment (peak_plot_EVT, d, org.Hs.eg.db, enrich_area='KEGG', log_FC_thres=thres)
display_cluster_enrichment (kk_EVT, show_graph='emap', feature_vec=
                                     peak_plot_EVT$group, show_num=30) + labs (fill='stage')
#
#PI3K <- gene_per_term (kk_EVT, 'PI3K-Akt', org.Hs.eg.db, return_val=T)[[1]]
#PI3K_sel <- peak_plot_EVT [PI3K,] %>% slice_max(val, n=4) %>% dplyr::select (feature) %>% deframe ()
PI3K_sel <- c('FLT4', 'LAMA4', 'CSF1R', 'EFNA1')
#hif <- gene_per_term (kk_EVT, 'HIF-1', org.Hs.eg.db, return_val=T)[[1]]
#hif_sel <- peak_plot_EVT [hif,] %>% slice_max(val, n=4) %>% dplyr::select (feature) %>% deframe ()
hif_sel <- c('TIMP1', 'TFRC', 'VEGFA', 'EGFR')

# obtain expresion matrix
p5 <- gene_over_pseudotime (pred_all, exp_mat, JAK_sel, show_meta,
                            color_feature = 'broad_type', num_col=2,
                            )+ggtitle ('JAK-STAT')
p6 <- gene_over_pseudotime (pred_all, exp_mat, PI3K_sel, show_meta,
                            color_feature = 'broad_type', num_col=2
                            )+ggtitle ('PI3K-AKT')
p7 <- gene_over_pseudotime (pred_all, exp_mat, hif_sel, show_meta,
                            color_feature = 'broad_type', num_col=2,
                            )+ggtitle ('HIF-1')
p567 <- lapply (list (p5, p6, p7), function (gx){gx + coord_cartesian (clip='off')+labs(color='cell type')})
p567 <- ggpubr::ggarrange (plotlist=p567, nrow=1, common.legend=T, legend='right')

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
sel_genes <- rownames (vivo_scale)
names (sel_genes) <- c('TB', 'non-specific', 'STB', 'non-specific', 'EPI', 'CTB', 'ICM', 'ICM',
                       'non-specific', 'EVT', 'cleavage')
row_levels <- partial_relevel (names (sel_genes)) %>% levels()
p8 <- seurat_heat (vivo_scale, group.by=c('broad_type', 'date'),
                 row_scale=T, color_row= sel_genes,
                 column_reorder_levels = list (format_conf$branch_order, format_conf$cell_order),
                 column_legend_labels= c('cell type', 'date'), 
                 row_legend_labels='WGCNA clusters',
                 cluster_rows=T, heat_name='module score', center_scale=T,
                 column_rotation=90,
                 main_width=16, main_height=11,
                 automatic=F, left_HA=F, slot_data='counts'
)

# ----------figure I: ICM specific module----------
save_dir4 <- paste (root, 'manuscript/figure4', sep='/')
color_row <- read.csv ( paste ( save_dir4,  'WGCNA/module_genes.csv' , sep='/'), row.names=1)
gene_list <- lapply (as.list (colnames (color_row) ), function (x) {
                             unique (color_row [, x]) })
names (gene_list) <- colnames (color_row)

# load data
data (TF) #load a vector of TF gene names
TF_WG <- all_data[TF,all_data$date != 'in_vitro']
fil_vivo <- filter_genes (TF_WG, 0.2)
rm (TF_WG)

save_dir1 <- paste (root, 'manuscript/figure1', sep='/')
markers <- find_DE_genes (all_data, save_dir1, group.by='broad_type', label='all_vivo')
merged7_8 <- c(as.character (gene_list$GC7), as.character (gene_list$GC8))
plot_gene <- list (ICM=merged7_8, TB=gene_list$GC1, STB=gene_list$GC3, EVT=gene_list$GC10)
plotlist <- custom_net_diff_nets (fil_vivo, plot_gene, markers, nudge_ratio=0.1, size_thres=0.2)
p9 <- plotlist[[1]]

# combine all the plots
p34_list <- list (theme (legend.position='right', legend.box='vertical'),
                  labs(color='cell type'))
grob_list <- list (p1+theme (legend.position='top'), 
                   p2+theme (legend.position='top'), 
                   p3+p34_list, 
                   p4+p34_list, 
                   p567, p8, p9)
lay_mat <- matrix(c(1, 1, 3, 3, 4, 4, 
                    2, 2, 3, 3, 4, 4,
                    5, 5, 5, 5, 5, 5,
                    6, 6, 6, 6, 7, 7
                    ),
                  nrow=6) %>% t()
arrange_plots (grob_list, paste (sup_save_dir, 'final_figureS2.pdf', sep='/'), 
                  lay_mat, plot_width=3, plot_height=6.3, margin_width=1)

save_indiv_plots (grob_list, paste (sup_save_dir, 'figureS2', sep='/'),
                  lay_mat, plot_width=3, plot_height=7
)
