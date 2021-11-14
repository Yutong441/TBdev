# WGCNA over pseudotime
setwd ('..')
devtools::load_all()
library (tidyverse)
library (org.Hs.eg.db)
library (GOSemSim)

root_dir <- '/mnt/c/Users/Yutong/Documents/bioinformatics/reproduction/'
root <- paste (root_dir, 'results/', sep='/')
merge_dir <- paste (root, 'XLYBPZ_Dylan_dir', sep='/')
save_dir <- paste (root, 'manuscript/figure3', sep='/')
sup_save_dir <- paste (root, 'manuscript/figureS3', sep='/')
x <- load (paste (merge_dir, 'final_merged_tb.Robj', sep='/') )
all_data <- get (x)

# ----------figure S3A-B: 2-WD ----------
data_dir <- paste (root_dir, 'GPLVM/hier_GP_tf1/', sep='/')
pred_all <- read.csv (paste (data_dir, 'result/prediction_matern_500.csv', sep='/' ))
pred_all$branch <- c('EVT_b', 'STB_b')[as.factor (pred_all$branch)]
exp_mat <- read.csv (paste (data_dir, 'data/STREAM_data.csv', sep='/' ), row.names=1)
pseudotime <- read.csv (paste (data_dir, 'result/infer_pt_matern.csv', sep='/' ), row.names=1)
exp_mat %>% as.matrix () %>% t () %>% scale () %>% data.frame () %>% 
        add_column (pseudotime=pseudotime [,'pt_mean']) -> exp_mat

markers <- get_DE_from_KL (pred_all, 'EVT_b', 'STB_b', divergence='WD')
markers$sign_diver <- sign (markers$logFC)*markers$divergence

show_meta <- all_data@meta.data [!is.na (all_data$MGP_PT) & !all_data$broad_type %in% c('EPI', 'PE') ,]
markers %>% slice_max (sign_diver, n=9) %>% dplyr::select (feature) %>% deframe() -> DE_genes1
p1 <- gene_over_pseudotime (pred_all, exp_mat, DE_genes1, show_meta, color_feature = 'broad_type', num_col=3)
markers %>% slice_min (sign_diver, n=9) %>% dplyr::select (feature) %>% deframe() -> DE_genes2
p2 <- gene_over_pseudotime (pred_all, exp_mat, DE_genes2, show_meta, color_feature = 'broad_type', num_col=3)

# ----------figure S3C----------
markers$group [markers$sign_diver <0 ] <- 'STB_b'
markers$group [markers$sign_diver ==0 ] <- 'none'
markers$logFC <- markers$divergence
thres <- quantile (markers$logFC, 0.7)

d <- godata(org.Hs.eg.db, ont="BP")
kk <- compare_cluster_enrichment (markers, d, org.Hs.eg.db, enrich_area='KEGG', log_FC_thres=thres)
set.seed (100)
p3 <- display_cluster_enrichment (kk, show_graph='emap', feature_vec=
                                     markers$group, show_num=20) + labs (fill='')

# ----------Module score over pseudotime----------
# figure S3E
sup_save_dir2 <- paste (root, 'manuscript/figureS2', sep='/')
color_row <- read.csv ( paste ( sup_save_dir2, 'WGCNA/module_genes.csv' , sep='/'), row.names=1)
gene_list <- lapply (as.list (colnames (color_row) ), function (x) {unique (color_row [, x]) })
names (gene_list) <- colors2labels (colnames (color_row), prefix='GC')

module_score <- get_module_score (all_data, paste (sup_save_dir2, 
                'WGCNA/Data_module_score.csv', sep='/'), pgenes=gene_list)
colnames (module_score) <- gsub ('^X', '', colnames (module_score))
module_score <- module_score [, match (rownames (show_meta), colnames (module_score) ) ]
module_data <- Seurat::CreateSeuratObject ( module_score, meta.data = show_meta)
module_data <- module_data [, !is.na (module_data$MGP_PT) ]
sel_genes <- rownames (module_data)
names (sel_genes) <- c('ECM', 'steroid', 'metabolism', 'ECM', 'NA', 'cytokine',
                       'immune', 'NA', 'immune', 'pluripotency',
                       'pluripotency')

p5 <- seurat_heat (module_data, group.by=c('epil_branch','broad_type'),
                 row_scale=T, color_row= sel_genes,
                 column_reorder_levels = list (CT$branch_order, CT$cell_order),
                 column_legend_labels= c('branch', 'cell type'), 
                 row_legend_labels='WGCNA clusters',
                 cluster_rows=T, heat_name='norm count', center_scale=T,
                 group_order = order (module_data$MGP_PT) )

# figure S3D: GO/KEGG on WGCNA
kk_wg <- compare_cluster_enrichment (gene_list, d, org.Hs.eg.db, enrich_area='KEGG')
color_cluster <- colnames (color_row)
names (color_cluster) <- names(gene_list)
set.seed (100)
p4 <- display_cluster_enrichment (kk_wg, show_graph='emap', feature_vec=names (gene_list), 
                                     show_num=20) + labs (fill = '') +
                                scale_fill_manual (values=color_cluster)
# arrange figures
grob_list <- list (p1, p2, p3, p4, p5)
lay_mat <- cbind (c(1,3,5), c(2,4,5))
arrange_plots (grob_list, paste (sup_save_dir, 'final_figureS3.pdf', sep='/'),  lay_mat)

# find diverging genes
stb_mark <- markers %>% filter (group == 'STB_branch')
GP$save_peak_genes (stb_mark, paste (sup_save_dir2, 'PT_marker/WD_STB_MAPK.csv', sep='/'), 'MAPK', kk, max_num=50, arrange_val='divergence')
GP$save_peak_genes (stb_mark, paste (sup_save_dir2, 'PT_marker/WD_STB_PI3K.csv', sep='/'), 'PI3K', kk, max_num=50, arrange_val='divergence')
GP$save_peak_genes (stb_mark, paste (sup_save_dir2, 'PT_marker/WD_STB_all.csv', sep='/'), arrange_val='divergence')

evt_mark <- markers %>% filter (group == 'EVT_branch')
GP$save_peak_genes (evt_mark, paste (sup_save_dir2, 'PT_marker/WD_EVT_AGE.csv', sep='/'), 'AGE-RAGE', kk, max_num=50, arrange_val='divergence')
GP$save_peak_genes (evt_mark, paste (sup_save_dir2, 'PT_marker/WD_EVT_all.csv', sep='/'), arrange_val='divergence')
