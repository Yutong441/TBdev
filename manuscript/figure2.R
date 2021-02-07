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
new_name <- c('TB_stem', 'EVT_branch', 'STB_branch', NA)
old_name <- c('S2,S1', 'S0,S1', 'S3,S1', 'blank')
data (format_conf)
epg$branch_name <- partial_relevel (new_name[match (epg$branch, old_name )], format_conf$branch_order)

# obtain label data
epg %>% arrange (branch_name) %>% mutate( index = 1:nrow (epg)) %>% 
        group_by (branch) %>% slice_min (index, n=1) %>% drop_na() %>% 
        as.data.frame () %>% dplyr::select (x, y, z, branch_name) %>% 
        magrittr::set_colnames (c('x', 'y', 'z', 'feature')) -> label_epg

metadata %>% slice_min (PT3, n=nrow(metadata)-3) -> plot_met
p1 <- dim_red_3D_traj (plot_met, 'PT1', 'PT2', 'PT3', 'broad_type', epg, 'x', 'y',
                    'z', 'branch_name', all_theta=50, all_phi=0, further_repel=T,
                    repel_force=0.5, lab_just=c(0.08, 0.02, 0.02), magnify_text=1.3, 
                    label_traj_text=label_epg, hor_just=0.1, dim_vjust=4) + 
labs (fill='cell type')
p1

# ----------Figure 2C: temporal clustering for STB----------
data_dir <- paste (root_dir, 'GPLVM/hier_GP_tf1/', sep='/')
pred_all <- read.csv (paste (data_dir, 'result/prediction_matern_500.csv', sep='/' ))
pred_all$x <- pred_all$x - min (all_data$ori_MGP_PT, na.rm=T)
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
new_clust <- c ('stage I', 'stage II', rep ('stage III', 4), 'stage IV'  )
peak_plot$group <- factor (new_clust [match (peak_plot$cluster, old_clust) ],levels=unique (new_clust))
write.csv (peak_plot, paste (save_dir, 'Switch_STB.csv', sep='/'))

metadata <- all_data@meta.data [!is.na (all_data@meta.data$MGP_PT) & !all_data@meta.data$broad_type %in% c('EPI', 'PE'), ]
metadata_STB <- metadata [all_data$epil_branch != 'EVT_branch' ,]
metadata_STB <- metadata_STB [metadata_STB$broad_type != 'EVT',]

show_data <- all_data [, match (rownames (metadata), colnames (all_data) ) ]
plot_heat1 <- show_data [, show_data$epil_branch != 'EVT_branch']
devtools::load_all ('..', export_all = F)
sel_genes_stb <- peak_gene_for_heatmap (peak_plot, min_genes=2)

p2 <- seurat_heat (plot_heat1, sel_genes_stb, group.by=c('broad_type'),
                    row_scale=T, show_row_names=F, 
                    group_order=order(plot_heat1$MGP_PT),
                    column_legend_labels = 'cell type', column_split=NA,
                    row_legend_labels ='stage',heat_name = 'norm count',
                    row_reorder_levels = c('stage I', 'stage II', 'stage III', 'stage IV'),
                    automatic=F, center_scale=T,
                    main_width=2.5, main_height=14,
)

# ----------Figure 2D: temporal clustering for EVT----------
peak_plot_EVT <- find_transition_from_seurat (PT_list_inf[[1]])
set.seed (100)
peak_plot_EVT$peak_time %>% Mclust(G=1:20) -> d_clust
peak_plot_EVT$cluster<- paste ('cluster', d_clust$classification, sep='')
old_clust <- paste ('cluster', 1:7, sep='')
new_clust <- c (rep('stage I',2), rep('stage II',2), 'stage IV', 'stage III', 'stage IV')
peak_plot_EVT$group <- factor (new_clust [match (peak_plot_EVT$cluster, old_clust) ],levels=unique (new_clust))
write.csv (peak_plot_EVT, paste (save_dir, 'Switch_EVT.csv', sep='/'))

plot_heat2 <- show_data [, show_data$epil_branch != 'STB_branch' & show_data$broad_type != 'STB']
sel_genes_evt <- peak_gene_for_heatmap (peak_plot_EVT, min_genes=2)

p3 <- seurat_heat (plot_heat2, sel_genes_evt, group.by=c('broad_type'),
                    row_scale=T, show_row_names=F, 
                    group_order=order(plot_heat2$MGP_PT),
                    column_legend_labels = 'cell type', column_split=NA,
                    row_legend_labels ='stage',heat_name = 'norm count',
                    row_reorder_levels = c('stage I', 'stage II', 'stage III', 'stage IV'),
                    automatic=F, center_scale=T,
                    main_width=2.5, main_height=14,
)

# ----------Figure 2E-F: KEGG----------
peak_plot$logFC <- peak_plot$val
thres <- quantile (peak_plot$logFC, 0.7)
kk_stb <- compare_cluster_enrichment (peak_plot, d, org.Hs.eg.db,
                                      enrich_area='KEGG', log_FC_thres=thres)
p4 <- display_cluster_enrichment (kk_stb, show_graph='emap', feature_vec=
                                     peak_plot$group, show_num=20) + labs (fill='')

peak_plot_EVT$logFC <- peak_plot_EVT$val
thres <- quantile (peak_plot_EVT$logFC, 0.7)
kk_evt <- compare_cluster_enrichment (peak_plot_EVT, d, org.Hs.eg.db,
                                      enrich_area='KEGG', log_FC_thres=thres)
p5 <- display_cluster_enrichment (kk_evt, show_graph='emap', feature_vec=
                                     peak_plot_EVT$group, show_num=20) + labs (fill='')
# ----------Figure 2G: pathway module score along pseudotime----------
module_scores <- get_module_score (all_data,
                              save_path=paste (sup_save_dir, 'Data_module_scores.csv', sep='/'))

module_scores <- module_scores [, match (rownames (metadata), colnames (module_scores) ) ]
rownames (module_scores) <- gsub ('\\.', '-', rownames (module_scores) )
path_seurat <- Seurat::CreateSeuratObject (module_scores, meta.data=metadata)
data (format_conf)
p6 <- seurat_heat (path_seurat, group.by=c('epil_branch','broad_type'),
                 color_row=rownames (path_seurat), row_scale=T,
                 column_reorder_levels = list (format_conf$branch_order, format_conf$cell_order),
                 column_legend_labels= c('branch', 'cell type'), 
                 heat_name = 'norm count', center_scale=T,
                 group_order = order (path_seurat$MGP_PT), automatic=F,
                 show_column_bars = c(F,T), cluster_rows=T
)

# ----------Figure 2B: schematic----------
# metadata for plotting the genes
save_dir <- paste (root, 'manuscript/figure2', sep='/')
peak_EVT <- read.csv (paste (save_dir, 'Switch_EVT.csv', sep='/'))
peak_STB <- read.csv (paste (save_dir, 'Switch_STB.csv', sep='/'))
peak_EVT$epil_branch <- 'EVT_branch'
peak_STB$epil_branch <- 'STB_branch'
peak_df <- rbind (peak_EVT, peak_STB)
pg <- graph_abs_time (metadata, peak_df, save_dir=paste (save_dir, 'graph_abs', sep='/'))

# ----------Figure H-K: WGCNA nets----------
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
plot_gene <- list (ICM=merged7_8, TB=gene_list$GC1, CTB=gene_list$GC6, STB=gene_list$GC3, EVT=gene_list$GC10)

plotlist <- custom_net_diff_nets (fil_vivo, plot_gene, markers, nudge_ratio=0.3, size_thres=0.2)
plotlist2 <- lapply (plotlist, function(gx){gx+theme(aspect.ratio=0.85)})
plotlist2 [[2]] <- plotlist[[2]] + ggtitle ('pre-implant')

# arrange figures
grob_list <- c(list (p1+theme (aspect.ratio=0.7), 
                     pg, p2, p3, p4+theme(aspect.ratio=0.85), 
                     p5+theme(aspect.ratio=0.85), 
                     p6
                     ), plotlist2)

lay_mat <- matrix(c(1, 1, 1, 1, 2, 2, 
                    1, 1, 1, 1, 2, 2,
                    3, 3, 4, 4, 5, 5,
                    3, 3, 4, 4, 6, 6,
                    7, 7, 7, 7, 8, 8,
                    7, 7, 7, 7, 9, 9,
                    10, 10, 11, 11, 12,12 
                    ),
                  nrow=6) %>% t()
arrange_plots (grob_list, paste (save_dir, 'final_figure2.pdf', sep='/'), 
                  lay_mat, plot_width=3, plot_height=3.5)
save_indiv_plots (grob_list, paste (save_dir, 'figure2', sep='/'), 
                  lay_mat, plot_width=3, plot_height=3.5)

# save key genes
#save_peak_genes (peak_plot, paste (sup_save_dir, 'PT_marker/Peak_STB_MAPK.csv', sep='/'), 'MAPK', kk, max_num=50)
#save_peak_genes (peak_plot, paste (sup_save_dir, 'PT_marker/Peak_STB_JAK.csv', sep='/'), 'JAK-STAT', kk, max_num=50)
#save_peak_genes (peak_plot, paste (sup_save_dir, 'PT_marker/Peak_STB_all.csv', sep='/'))
