devtools::load_all('..', export_all=F)
library (tidyverse)
library (GOSemSim)
library (org.Hs.eg.db)

# load data
root_dir <- '/mnt/c/Users/Yutong/Documents/bioinformatics/reproduction/'
root <- paste (root_dir, 'results/', sep='/')
merge_dir <- paste (root, 'XLYBPZ_Dylan_dir', sep='/')
sup_save_dir <- paste (root, 'manuscript/figureS2', sep='/')
all_data <- get (load (paste (merge_dir, 'final_merged_tb.Robj', sep='/') ))

# ----------figure S2A----------
show_meta <- all_data@meta.data [!(is.na(all_data$MGP_PT) | all_data$broad_type %in% c('PE', 'EPI')), ]
p1 <- pseudo_real_time (show_meta, 'date', 'MGP_PT', 'broad_type')

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
p4 <- gene_over_pseudotime (pred_all, exp_mat, DE_genes1, show_meta, color_feature = 'broad_type', num_col=3)
markers %>% slice_min (sign_diver, n=9) %>% dplyr::select (feature) %>% deframe() -> DE_genes2
p3 <- gene_over_pseudotime (pred_all, exp_mat, DE_genes2, show_meta, color_feature = 'broad_type', num_col=3)

# ----------figure S3B----------
markers$group [markers$sign_diver <0 ] <- 'STB_branch'
markers$group [markers$sign_diver ==0 ] <- 'none'
markers$group <- partial_relevel (markers$group)
markers$logFC <- markers$divergence
thres <- quantile (markers$logFC, 0.7)

d <- godata(org.Hs.eg.db, ont="BP")
kk <- compare_cluster_enrichment (markers, d, org.Hs.eg.db, enrich_area='KEGG', log_FC_thres=thres)
p2 <- display_cluster_enrichment (kk, show_graph='emap', feature_vec=
                                     markers$group, show_num=20) + labs (fill='')
p2

# ----------figure S2E----------
save_dir4 <- paste (root, 'manuscript/figure4', sep='/')
color_row <- read.csv ( paste ( save_dir4, 
                        'WGCNA/module_genes.csv' , sep='/'), row.names=1)
gene_list <- lapply (as.list (colnames (color_row) ), function (x) {
                             unique (color_row [, x]) })
ori_names <- colnames (color_row)
names (gene_list) <- colors2labels (colnames (color_row), prefix='GC')

module_score <- get_module_score (invivo, paste (save_dir4, 
                'WGCNA/Data_module_score.csv', sep='/'), pgenes=gene_list)
colnames (module_score) <- gsub ('^X', '', colnames (module_score))

show_meta <- all_data@meta.data [!is.na (all_data$MGP_PT) & !all_data$broad_type %in% c('EPI', 'PE') ,]
module_score <- module_score [, match (rownames (show_meta), colnames (module_score) ) ]
module_data <- Seurat::CreateSeuratObject ( module_score, meta.data = show_meta)
module_data <- module_data [, !is.na (module_data$MGP_PT) ]
sel_genes <- rownames (module_data)
names (sel_genes) <- c('TB', 'EVT', 'STB', 'CTB', 'non-specific', 
                       'CTB', 'cleavage', 'cleavage', 'non-specific', 'EVT',
                       'non-specific')
row_levels <- partial_relevel (names (sel_genes)) %>% levels()

data (format_conf)
devtools::load_all('..', export_all=F)
p5 <- seurat_heat (module_data, group.by=c('epil_branch','broad_type'),
                 row_scale=T, color_row= sel_genes,
                 column_reorder_levels = list (format_conf$branch_order, format_conf$cell_order),
                 column_legend_labels= c('branch', 'cell type'), 
                 show_column_bars=c(F,T),
                 row_legend_labels='WGCNA clusters',
                 cluster_rows=T, heat_name='norm count', center_scale=T,
                 group_order = order (module_data$MGP_PT),
                 automatic=F, left_HA=F
)

grob_list <- list (p1, p2, p3, p4, p5)
lay_mat <- matrix(c(1, 2, 
                    3, 4,
                    5, 5
                    ),
                  nrow=2) %>% t()
arrange_plots (grob_list, paste (sup_save_dir, 'final_figureS2.pdf', sep='/'), 
                  lay_mat, plot_width=7, plot_height=7)

# find peak genes
save_peak_genes (peak_plot, paste (sup_save_dir, 'PT_marker/Peak_EVT_PI3K.csv', sep='/'), 'PI3K-Akt', kk, max_num=50)
save_peak_genes (peak_plot, paste (sup_save_dir, 'PT_marker/Peak_EVT_HIF1.csv', sep='/'), 'HIF-1', kk, max_num=50)
save_peak_genes (peak_plot, paste (sup_save_dir, 'PT_marker/Peak_EVT_relaxin.csv', sep='/'), 'Relaxin', kk, max_num=50)
save_peak_genes (peak_plot, paste (sup_save_dir, 'PT_marker/Peak_EVT_all.csv', sep='/'))
