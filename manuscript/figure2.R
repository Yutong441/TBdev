setwd ('..')
library (tidyverse)
library (mclust)
library (GOSemSim)
devtools::load_all()

root_dir <- '/mnt/c/Users/Yutong/Documents/bioinformatics/reproduction/'
root <- paste (root_dir, 'results/', sep='/')
merge_dir <- paste (root, 'XLYBPZ_Dylan_dir', sep='/')
save_dir <- paste (root, 'manuscript/figure2', sep='/')
sup_save_dir <- paste (root, 'manuscript/figureS2', sep='/')
x <- load (paste (merge_dir, 'final_merged_tb.Robj', sep='/') )
all_data <- get (x)

# ----------Figure 2A: Epilgraph----------
epg <- read.csv (paste (sup_save_dir, 'result/STREAM_graph.csv', sep='/'))
metadata <- all_data@meta.data 
metadata %>% filter (!broad_type %in% c('EPI', 'PE')) %>% filter (!is.na (PT1)) -> metadata
p1 <- dim_red_3D_traj (metadata, 'PT1', 'PT2', 'PT3', 'broad_type', epg, 'x', 'y',
                    'z', 'branch', all_theta=50, all_phi=0, further_repel=T,
                    repel_force=0.5, lab_just=c(0.08, 0.02, 0.02)) + labs (fill='')
# ----------Figure 2B: temporal clustering----------
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
metadata_STB <- metadata [all_data$epil_branch != 'EVT_b' ,]
metadata_STB <- metadata_STB [metadata_STB$broad_type != 'EVT',]
set.seed (100)
p2 <- time_cluster_plot (peak_plot, metadata_STB, color_by='group', show_text_prop=0.999, 
                         repel_force=10, repel_point=0.98) 

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
kk <- compare_cluster_enrichment (peak_plot, d, enrich_area='KEGG', log_FC_thres=thres)
display_cluster_enrichment (kk, show_graph='emap', feature_vec=
                                     peak_plot$group, show_num=30) + labs (fill='stage')
one_term <- gene_per_term (kk, 'MAPK', return_val=T)

# obtain expresion matrix
exp_mat <- read.csv (paste (data_dir, 'data/STREAM_data.csv', sep='/' ), row.names=1)
pseudotime <- read.csv (paste (data_dir, 'result/infer_pt_matern.csv', sep='/' ), row.names=1)
exp_mat %>% as.matrix () %>% t () %>% scale () %>% data.frame () %>% 
        add_column (pseudotime=pseudotime [,'pt_mean']) -> exp_mat

pred_all$branch <- c('EVT_b', 'STB_b')[as.factor (pred_all$branch)]
p5 <- gene_over_pseudotime (pred_all, exp_mat, one_term[[1]], metadata,
                            color_feature = 'broad_type', num_col=9,
                            peak_data=peak_plot)

# arrange figures
grob_list <- list (p1, p2+labs (fill='cell type', color='stage'), 
                   p3, p4+ theme (legend.position='right',aspect.ratio=0.5), p5)
lay_mat <- matrix(c(1, 1, 2, 2, 2, 2, 
                    3, 3, 3, 3, 4, 4,
                    5, 5, 5, 5, 5, 5),
                  nrow=6) %>% t()
arrange_plots (grob_list, paste (save_dir, 'final_figure2.pdf', sep='/'), 
                  lay_mat, plot_width=3., plot_height=6)

# save key genes
#save_peak_genes (peak_plot, paste (sup_save_dir, 'PT_marker/Peak_STB_MAPK.csv', sep='/'), 'MAPK', kk, max_num=50)
#save_peak_genes (peak_plot, paste (sup_save_dir, 'PT_marker/Peak_STB_JAK.csv', sep='/'), 'JAK-STAT', kk, max_num=50)
#save_peak_genes (peak_plot, paste (sup_save_dir, 'PT_marker/Peak_STB_all.csv', sep='/'))
