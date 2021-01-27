#library (TBdev)
devtools::load_all ('..', export_all = F)
library (tidyverse)
library (mclust)
library (GOSemSim)
library (org.Hs.eg.db)

root_dir <- '/mnt/c/Users/Yutong/Documents/bioinformatics/reproduction/'
root <- paste (root_dir, 'results/', sep='/')
merge_dir <- paste (root, 'XLYBPZ_Dylan_dir', sep='/')
save_dir <- paste (root, 'manuscript/figure2', sep='/')
sup_save_dir <- paste (root, 'manuscript/figureS2', sep='/')
all_data <- get (load (paste (merge_dir, 'final_merged_tb.Robj', sep='/') ))

# ----------Figure 2A: Epilgraph----------
epg <- read.csv (paste (sup_save_dir, 'result/STREAM_graph.csv', sep='/'))
metadata <- all_data@meta.data 
metadata %>% filter (!broad_type %in% c('EPI', 'PE')) %>% filter (!is.na (PT1)) -> metadata
new_name <- c('main', 'EVT_branch', 'STB_branch', NA)
old_name <- c('S2,S1', 'S0,S1', 'S3,S1', 'blank')
data (format_conf)
epg$branch_name <- partial_relevel (new_name[match (epg$branch, old_name )], format_conf$branch_order)

# obtain label data
epg %>% arrange (branch_name) %>% mutate( index = 1:nrow (epg)) %>% 
        group_by (branch) %>% slice_min (index, n=1) %>% drop_na() %>% 
        as.data.frame () %>% dplyr::select (x, y, z, branch_name) %>% 
        magrittr::set_colnames (c('x', 'y', 'z', 'feature')) -> label_epg

p1 <- dim_red_3D_traj (metadata, 'PT1', 'PT2', 'PT3', 'broad_type', epg, 'x', 'y',
                    'z', 'branch_name', all_theta=50, all_phi=0, further_repel=T,
                    repel_force=0.5, lab_just=c(0.08, 0.02, 0.02), magnify_text=1.3, label_traj_text=label_epg) + 
labs (fill='cell type')

# ----------Figure 2B: temporal clustering for STB----------
data_dir <- paste (root_dir, 'GPLVM/hier_GP_tf1/', sep='/')
pred_all <- read.csv (paste (data_dir, 'result/prediction_matern_500.csv', sep='/' ))
label_list <- list ( c('branch0', 'branch1'),  c('branch0', 'branch2') )
PT_list_inf <- raw_to_seurat (pred_all, label_list)
peak_plot <- find_transition_from_seurat (PT_list_inf[[2]])

# It is difficult to automate this function because I do not know the
# clustering results a priori
# Set seed is very important here
set.seed (100)
peak_plot$peak_time  %>% Mclust(G=1:20) -> d_clust
peak_plot$cluster<- paste ('cluster', d_clust$classification, sep='')
old_clust <- paste ('cluster', 1:7, sep='')
new_clust <- c ('early', 'intermediate1', rep ('intermediate2', 4), 'advanced'  )
peak_plot$group <- factor (new_clust [match (peak_plot$cluster, old_clust) ],levels=unique (new_clust))

metadata <- all_data@meta.data [!is.na (all_data@meta.data$MGP_PT) & !all_data@meta.data$broad_type %in% c('EPI', 'PE'), ]
metadata_STB <- metadata [all_data$epil_branch != 'EVT_branch' ,]
metadata_STB <- metadata_STB [metadata_STB$broad_type != 'EVT',]
p2 <- time_cluster_plot (peak_plot, metadata_STB, color_by='group', show_text_prop=0.999, 
                         repel_force=10, repel_point=0.98) 
# ----------Figure 2B: temporal clustering for EVT----------
peak_plot_EVT <- find_transition_from_seurat (PT_list_inf[[1]])
set.seed (100)
peak_plot_EVT$peak_time  %>% Mclust(G=1:20) -> d_clust
peak_plot_EVT$cluster<- paste ('cluster', d_clust$classification, sep='')
old_clust <- paste ('cluster', 1:7, sep='')
new_clust <- c (rep('early',2), rep('intermediate1',2), 'advanced', 'intermediate2', 'advanced')
peak_plot_EVT$group <- factor (new_clust [match (peak_plot_EVT$cluster, old_clust) ],levels=unique (new_clust))

metadata_EVT <- metadata [metadata $epil_branch != 'STB_branch',]
metadata_EVT <- metadata_EVT [metadata_EVT$broad_type != 'STB',]

# ----------Figure 2C: pathway module score along pseudotime----------
module_scores <- get_module_score (all_data,
                              save_path=paste (sup_save_dir, 'Data_module_scores.csv', sep='/'))

module_scores <- module_scores [, match (rownames (metadata), colnames (module_scores) ) ]
rownames (module_scores) <- gsub ('\\.', '-', rownames (module_scores) )
path_seurat <- Seurat::CreateSeuratObject (module_scores, meta.data=metadata)
data (format_conf)
p3 <- seurat_heat (path_seurat, group.by=c('epil_branch','broad_type'),
                 color_row=rownames (path_seurat), row_scale=T,
                 column_reorder_levels = list (format_conf$branch_order, format_conf$cell_order),
                 column_legend_labels= c('branch', 'cell type'), 
                 heat_name = 'norm count', center_scale=T,
                 group_order = order (path_seurat$MGP_PT), automatic=F)

# ----------Figure 2D: cell cycle over time----------
exp_mat_cc <- as.matrix (Seurat::GetAssayData (all_data, assay='RNA', slot='data'))
ans <- get_phase_score (exp_mat_cc)
ans <- data.frame ( ans [ match( rownames (metadata), rownames (ans)  ), ] )
plot_data <- cbind (ans, metadata)
p4 <- cycle_over_time (plot_data, 'broad_type', time_col='MGP_PT') + xlab ('pseudotime')

# ----------Figure 2E: MAPK over time----------
d <- godata(org.Hs.eg.db, ont="BP")
peak_plot$logFC <- peak_plot$val
thres <- quantile (peak_plot$logFC, 0.7)
kk <- compare_cluster_enrichment (peak_plot, d, org.Hs.eg.db, enrich_area='KEGG', log_FC_thres=thres)
display_cluster_enrichment (kk, show_graph='emap', feature_vec=
                                     peak_plot$group, show_num=30) + labs (fill='stage')
# genes in MAPK pathway
MAPK <- gene_per_term (kk, 'MAPK', org.Hs.eg.db, return_val=T)[[1]]
#peak_plot [MAPK,] %>% arrange (desc (val))
MAPK_sel <- c('GADD45G', 'RAP1B', 'DUSP8', 'RASA1')

# genes in JAK-STAT pathway
JAK <- gene_per_term (kk, 'JAK-STAT', org.Hs.eg.db, return_val=T)[[1]]
#peak_plot [JAK,] %>% arrange (desc (val))
JAK_sel <- peak_plot [JAK,] %>% slice_max(val, n=4) %>% dplyr::select (feature) %>% deframe ()

# PI3K-AKT
peak_plot_EVT$logFC <- peak_plot_EVT$val
thres <- quantile (peak_plot_EVT$logFC, 0.7)
kk_EVT <- compare_cluster_enrichment (peak_plot_EVT, d, org.Hs.eg.db, enrich_area='KEGG', log_FC_thres=thres)
display_cluster_enrichment (kk_EVT, show_graph='emap', feature_vec=
                                     peak_plot_EVT$group, show_num=30) + labs (fill='stage')

PI3K <- gene_per_term (kk_EVT, 'PI3K-Akt', org.Hs.eg.db, return_val=T)[[1]]
#peak_plot [PI3K,] %>% arrange (desc (val))
PI3K_sel <- peak_plot_EVT [PI3K,] %>% slice_max(val, n=4) %>% dplyr::select (feature) %>% deframe ()

# obtain expresion matrix
exp_mat <- read.csv (paste (data_dir, 'data/STREAM_data.csv', sep='/' ), row.names=1)
pseudotime <- read.csv (paste (data_dir, 'result/infer_pt_matern.csv', sep='/' ), row.names=1)
exp_mat %>% as.matrix () %>% t () %>% scale () %>% data.frame () %>% 
        add_column (pseudotime=pseudotime [,'pt_mean']) -> exp_mat

pred_all$branch <- c('EVT_branch', 'STB_branch')[as.factor (pred_all$branch)]

p5 <- gene_over_pseudotime (pred_all, exp_mat, MAPK_sel, metadata,
                            color_feature = 'broad_type', num_col=2,
                            peak_data=peak_plot)
p6 <- gene_over_pseudotime (pred_all, exp_mat, JAK_sel, metadata,
                            color_feature = 'broad_type', num_col=2,
                            peak_data=peak_plot)
p7 <- gene_over_pseudotime (pred_all, exp_mat, PI3K_sel, metadata,
                            color_feature = 'broad_type', num_col=2,
                            peak_data=peak_plot_EVT)

# arrange figures
grob_list <- list (p1, p2+labs (fill='cell type', color='stage')+theme(aspect.ratio=0.8), 
                   p3, p4+ theme (legend.position='right',aspect.ratio=0.5), p5,p6,p7)
lay_mat <- matrix(c(1, 1, 1, 2, 2, 2, 
                    3, 3, 3, 3, 4, 4,
                    5, 5, 6, 6, 7, 7),
                  nrow=6) %>% t()
arrange_plots (grob_list, paste (save_dir, 'final_figure2.pdf', sep='/'), 
                  lay_mat, plot_width=3., plot_height=6)

# save key genes
#save_peak_genes (peak_plot, paste (sup_save_dir, 'PT_marker/Peak_STB_MAPK.csv', sep='/'), 'MAPK', kk, max_num=50)
#save_peak_genes (peak_plot, paste (sup_save_dir, 'PT_marker/Peak_STB_JAK.csv', sep='/'), 'JAK-STAT', kk, max_num=50)
#save_peak_genes (peak_plot, paste (sup_save_dir, 'PT_marker/Peak_STB_all.csv', sep='/'))

# quantify switch point
peak_plot %>% filter (feature %in% one_term[[1]]) %>% summarise(pt = mean (peak_time)) 
# 1 -0.01059811
metadata_STB %>% group_by (broad_type) %>% summarise (min_pt = min(MGP_PT), max_pt = max(MGP_PT))

peak_EVT <- find_transition_from_seurat (PT_list_inf[[1]])
peak_EVT %>% filter (feature %in% one_term[[1]]) %>% summarise(pt = mean (peak_time)) 
# 1 -0.2877312

t.test (peak_EVT [match (one_term[[1]], peak_EVT$feature), 'peak_time'],
        peak_plot [match (one_term[[1]], peak_EVT$feature), 'peak_time'], paired=T
)

