setwd ('..')
devtools::load_all()
library (tidyverse)
library (GOSemSim)
library (org.Hs.eg.db)
library (mclust)

# load data
root_dir <- '/mnt/c/Users/Yutong/Documents/bioinformatics/reproduction/'
root <- paste (root_dir, 'results/', sep='/')
merge_dir <- paste (root, 'XLYBPZ_Dylan_dir', sep='/')
sup_save_dir <- paste (root, 'manuscript/figureS2', sep='/')
x <- load (paste (merge_dir, 'final_merged_tb.Robj', sep='/') )
all_data <- get (x)

# ----------figure S2A: transition markers----------
marker_dir <- paste (sup_save_dir, 'result/transition_markers/', sep='/') 
all_files <- list.files (marker_dir)
branch_name <- gsub ('^transition_markers_', '', all_files)
branch_name <- gsub ('.tsv$', '', branch_name)
markers <- list ()
for (i in 1:length (all_files)){
        one_data <- read.table (paste (marker_dir, all_files[i], sep='/') )
        colnames (one_data)[2] <- 'logFC'
        one_data$feature <- rownames (one_data)
        one_data$group <- branch_name [i]
        markers [[i]] <- one_data
}

markers <- do.call (rbind, markers)
# assign meaningful names to branches
ori_names <- c( 'S2_S1', 'S1_S0', 'S1_S3')
new_names <- c('main','EVT_b', 'STB_b')
markers$group <- new_names [match (markers$group, ori_names) ]
markers$group <- factor (markers$group, levels=new_names)

unique_DE_genes (markers %>% filter (stat > 0.) , 10) %>% 
        arrange (desc (stat) ) %>%
        dplyr::select (group, feature) %>% 
        deframe () -> DE_genes

data (CT)
tb_data <- all_data [, !(all_data$broad_type %in% CT$non_TB_lineage | all_data$date == 'in_vitro')]
order_heat <- gsub('^D', '', tb_data$date) %>% as.numeric () %>% order ()
p1 <- seurat_heat (tb_data, color_row=DE_genes, group.by = c('broad_type'), 
                       slot='data', heat_name='norm count',
                       column_legend_labels=c('cell type'), row_scale=T,
                       row_legend_labels ='stage', center_scale=T,
                       row_reorder_levels = c('main', 'EVT_b', 'STB_b'),
                       column_rotation=0, 
                       group_order=order_heat, automatic=F)
p1
#ComplexHeatmap::draw (p1, align_heatmap_legend='heatmap_center')
# ----------figure S2B----------
d <- godata(org.Hs.eg.db, ont="BP")
markers %>% filter (logFC > 0.25) -> sel_markers
kk_trans <- compare_cluster_enrichment (sel_markers, d, org.Hs.eg.db, enrich_area='KEGG')
set.seed (100)
p2 <- display_cluster_enrichment (kk_trans, show_graph='emap', feature_vec=
                                     sel_markers$group) + labs (fill='')
p2

# ----------figure S2C----------
show_meta <- all_data@meta.data [!(is.na(all_data$MGP_PT) | all_data$broad_type %in% c('PE', 'EPI')), ]
p3 <- pseudo_real_time (show_meta, 'date', 'MGP_PT', 'broad_type')

# ----------figure S2D----------
data_dir <- paste (root_dir, 'GPLVM/hier_GP_tf1/', sep='/')
pred_all <- read.csv (paste (data_dir, 'result/prediction_matern_500.csv', sep='/' ))
label_list <- list ( c('branch0', 'branch1'),  c('branch0', 'branch2') )
PT_list_inf <- raw_to_seurat (pred_all, label_list)
peak_plot <- find_transition_from_seurat (PT_list_inf[[1]])

set.seed (100)
peak_plot$peak_time  %>% Mclust(G=1:20) -> d_clust
peak_plot$cluster<- paste ('cluster', d_clust$classification, sep='')
old_clust <- paste ('cluster', 1:7, sep='')
new_clust <- c (rep('early',2), rep('intermediate1',2), 'advanced', 'intermediate2', 'advanced')
peak_plot$group <- factor (new_clust [match (peak_plot$cluster, old_clust) ],levels=unique (new_clust))

metadata_EVT <- show_meta [show_meta$epil_branch != 'STB_b',]
metadata_EVT <- metadata_EVT [metadata_EVT$broad_type != 'STB',]
p4 <- time_cluster_plot (peak_plot, metadata_EVT, color_by='group',
                         show_text_prop=0.998, repel_point=0.99, repel_force=1) 
p4
# ----------figure S2G: EVT temporal gene modules----------
peak_plot %>% dplyr::filter (val > quantile (peak_plot$val, 0.95) ) %>% 
        arrange (desc(val)) -> label_peak
label_peak %>% arrange (change_sign*val ) %>% dplyr::select (group, feature) %>% deframe () -> sel_genes
show_data <- all_data [, match (rownames (show_meta), colnames (all_data) ) ]
plot_heat1 <- show_data [, show_data$epil_branch != 'STB_b']
p7 <- seurat_heat (plot_heat1, sel_genes, group.by=c('broad_type'),
                    row_scale=T, show_row_names=F, 
                    group_order=order(plot_heat1$MGP_PT),
                    column_legend_labels = 'cell type', column_split=NA,
                    row_legend_labels ='stage',heat_name = 'norm count',
                    row_reorder_levels = c('early', 'intermediate1', 'intermediate2', 'advanced'),
                    automatic=F, center_scale=T,
)
p7

# ----------figure 2H----------
# figure 2E
peak_plot$logFC <- peak_plot$val
thres <- quantile (peak_plot$logFC, 0.7)
kk <- compare_cluster_enrichment (peak_plot, d, org.Hs.eg.db,
                                  enrich_area='KEGG', log_FC_thres=thres)
set.seed (100)
p8<- display_cluster_enrichment (kk, show_graph='emap', feature_vec=
                                     peak_plot$group, show_num=20) + labs (fill='')
path_val <- gene_per_term (kk, 'PI3K-Akt', org.Hs.eg.db, return_val=T)
path_val[[1]]
#  [1] "ANGPT2" "CCND1"  "COL4A1" "COL4A2" "CSF1R"  "CSH2"   "EFNA1"  "EGFR"   "FLT4"   "GH2"    "GNG11"  "GNG4"   "GNGT1"  "IL6"    "ITGA5"  "ITGB1" 
# [17] "JAK1"   "LAMA4"  "SPP1"   "VEGFA" 

# ----------figure S2E----------
peak_plot2 <- find_transition_from_seurat (PT_list_inf[[2]])
set.seed (100)
peak_plot2$peak_time  %>% Mclust(G=1:20) -> d_clust
peak_plot2$cluster<- paste ('cluster', d_clust$classification, sep='')

old_clust <- paste ('cluster', 1:7, sep='')
new_clust <- c ('early', 'intermediate1', rep ('intermediate2', 4), 'advanced'  )
peak_plot2$group <- new_clust [match (peak_plot2$cluster, old_clust) ]

peak_plot2 %>% dplyr::filter (val > quantile (peak_plot2$val, 0.95) ) %>% 
        arrange (desc(val)) -> label_peak
label_peak %>% arrange (change_sign*val ) %>% dplyr::select (group, feature) %>% deframe () -> sel_genes

plot_heat2 <- show_data [, show_data$epil_branch != 'EVT_b']
p5 <- seurat_heat (plot_heat2, sel_genes, group.by=c('broad_type'),
                    row_scale =T, show_row_names=F, 
                    group_order=order(plot_heat2$MGP_PT),
                    column_legend_labels = 'cell type', column_split=NA,
                    row_legend_labels ='stage', heat_name = 'norm count',
                    row_reorder_levels = c('early', 'intermediate1', 'intermediate2', 'advanced'),
                    automatic=F, center_scale=T,
)
p5

# ----------figure 2F----------
peak_plot2$logFC <- peak_plot2$val
thres <- quantile (peak_plot2$logFC, 0.7)
kk_stb <- compare_cluster_enrichment (peak_plot2, d, org.Hs.eg.db,
                                      enrich_area='KEGG', log_FC_thres=thres)
set.seed (100)
p6 <- display_cluster_enrichment (kk_stb, show_graph='emap', feature_vec=
                                     peak_plot$group, show_num=20) + labs (fill='')
path_val <- gene_per_term (kk_stb, 'MAPK', org.Hs.eg.db, return_val=T)
path_val[[1]]
#  [1] "ANGPT4"  "DDIT3"   "DUSP6"   "DUSP8"   "DUSP9"   "EPHA2"   "FGF18"   "GADD45B" "GADD45G" "GNG12"   "HSPA1B"  "HSPB1"   "IL1R1"   "IL1RAP"  "JUND"   
# [16] "MECOM"   "MKNK2"   "MYC"     "NR4A1"   "PDGFA"   "PGF"     "PPM1B"   "RAP1B"   "RASA1"   "RPS6KA5" "STMN1"   "TGFBR2" 

path_val <- gene_per_term (kk_stb, 'JAK-STAT', org.Hs.eg.db, return_val=T)
path_val[[1]]
#  [1] "CSF2RB" "CSH2"   "EGFR"   "GH2"    "IL2RB"  "JAK1"   "LEP"    "LIFR"   "OSMR"   "PRLR"  

# ----------arrange all----------
grob_list <- list (p1, p2, p3, p4+labs (color='stage', fill='cell type'), p5, p6, p7, p8)
lay_mat <- matrix(c(1, 1, 1, 1, 2, 2, 
                    3, 3, 4, 4, 4, 4,
                    5, 5, 5, 5, 6, 6,
                    7, 7, 7, 7, 8, 8),
                  nrow=6) %>% t()
arrange_plots (grob_list, paste (sup_save_dir, 'final_figureS2.pdf', sep='/'), 
                  lay_mat, plot_width=3., plot_height=6)

# find peak genes
save_peak_genes (peak_plot, paste (sup_save_dir, 'PT_marker/Peak_EVT_PI3K.csv', sep='/'), 'PI3K-Akt', kk, max_num=50)
save_peak_genes (peak_plot, paste (sup_save_dir, 'PT_marker/Peak_EVT_HIF1.csv', sep='/'), 'HIF-1', kk, max_num=50)
save_peak_genes (peak_plot, paste (sup_save_dir, 'PT_marker/Peak_EVT_relaxin.csv', sep='/'), 'Relaxin', kk, max_num=50)
save_peak_genes (peak_plot, paste (sup_save_dir, 'PT_marker/Peak_EVT_all.csv', sep='/'))
