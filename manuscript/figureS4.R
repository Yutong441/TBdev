# generate the figure S4 of the manuscript
# Identity mapping of Trophoblast stem cells

root_dir <- '/mnt/c/Users/Yutong/Documents/bioinformatics/reproduction/'
setwd (paste (root_dir, 'SingleCellR/utils', sep='/'))
library (modules)
library (tidyverse)
library (Seurat)
library (GOSemSim)
DIR <- modules::use ('dim_red.R')
ML <- modules::use ('marker_list.R')
DEG <- modules::use ('DE_gene.R')
PD <- modules::use ('plot_DR.R')
KG <- modules::use ('KEGG_path.R')
FO <- modules::use ('format.R')
WG <- modules::use ('WGCNA_utils.R')

root <- paste (root_dir, 'results/', sep='/')
merge_dir <- paste (root, 'XLYBPZ_Dylan_dir', sep='/')
save_dir <- paste (root, 'manuscript/figure3', sep='/')
sup_save_dir3 <- paste (root, 'manuscript/figureS3', sep='/')
sup_save_dir <- paste (root, 'manuscript/figureS4', sep='/')
x <- load (paste (merge_dir, 'final_merged_vivo.Robj', sep='/') )
all_data <- get (x)

# assign clusters
TSC_cells <- all_data$revised %in% c('hTSC_OKAE', 'hTSC_TURCO') 
all_data$new_cluster <- all_data$assigned_cluster

# ----------figure 3SA: diffusion map ----------
all_data2 <- all_data [, ! (all_data$broad_type %in% c( ML$non_emb_lineage,
                                                    ML$pre_imp_lineage))]

all_data2 <- DIR$run_dim_red (all_data2, run_diff_map=T, var_scale=T,
                             normalize=F, find_var_features=T, run_umap=F)
# filter out cells that interfere with aesthetics
show_data <- all_data2 [, !(all_data2$assigned_cluster == 'uCTB' ) ]
highlight <- show_data$new_cluster %in% ML$in_vitro_cells
p1 <- DIR$plot_dim_red (show_data, by_group = c('new_cluster'), DR='diff_map' , return_sep=T,
                    size_highlight=highlight, highlight_font=2, dims=c(1,2),
                    nudge_ratio=1.2, move_x=1.3, move_y=6)
p1[[1]]

# ----------figure S4B: over-representation for Okae----------
selected_cells <- c('eTB', 'iTB', 'aTB', 'vCTB1', 'hTSC_OKAE')
d <- godata('org.Hs.eg.db', ont="BP")
all_data$new_cluster2 <- as.character (all_data$assigned_cluster)
#all_data$new_cluster2 [TSC_cells] <- as.character (all_data$revised [TSC_cells])
OKAE_data <- all_data [, all_data$new_cluster2 %in% selected_cells]
OKAE_data$new_cluster2 <- ML$partial_relevel (OKAE_data$new_cluster2, ML$cell_order)
markers <- DEG$find_DE_genes (OKAE_data, sup_save_dir3, feature='new_cluster2', label='OKAE')
kk_okae <- KG$compare_cluster_enrichment (markers, d, enrich_area='KEGG')
p2<- KG$display_cluster_enrichment (kk_okae, show_graph='emap', feature_vec=
                                     markers$group, show_num=20) + labs (fill='cell type')
p2

# ----------figure S4C: for Turco----------
selected_cells <- c('vCTB1', 'vCTB2', 'vCTB3','hTSC_TURCO')
TURCO_data <- all_data [, all_data$new_cluster2 %in% selected_cells]
TURCO_data$new_cluster2 <- ML$partial_relevel (TURCO_data$new_cluster2, ML$cell_order)
markers <- DEG$find_DE_genes (TURCO_data, sup_save_dir3, feature='new_cluster2', label='TURCO')
kk_turco <- KG$compare_cluster_enrichment (markers, d, enrich_area='KEGG')
p3 <- KG$display_cluster_enrichment (kk_turco, show_graph='emap', feature_vec=
                                     markers$group, show_num=20) + labs (fill='cell type')

# ----------figure S4D: probability----------
sup_save_dir2 <- paste (root, 'manuscript/figureS2', sep='/')
epg <- read.csv (paste (sup_save_dir2, 'result/STREAM_graph.csv', sep='/'))
show_meta <- all_data@meta.data
show_meta %>% filter (!is.na (MGP_PT) ) %>% filter (!broad_type %in% c('PE', 'hTSC') ) -> show_meta

in_vitro <- c( 'hTSC_OKAE', 'hTSC_TURCO', 'hESC', 'hESC_YAN')
in_vivo <- as.character (unique (all_data$assigned_cluster) )
color_by <- unique (c(in_vivo, in_vitro))
p4 <- PD$dim_red_3D_traj (show_meta, 'PT1', 'PT2', 'PT3', color_by, epg, 'x',
                          'y', 'z', 'branch', all_theta=50, all_phi=0,
                          show_axes=T, show_label=F, num_col=9) +
                          labs (fill='probability')

# ----------figure S4E: pathway module----------
sup_save_dir2 <- paste (root, 'manuscript/figureS2', sep='/')
module_scores <- KG$get_module_score (all_data, all_path=KG$kg_pathway,
                                      save_path=paste (sup_save_dir2, 'Data_module_scores.csv', sep='/'))
colnames (module_scores) <- gsub ('^X', '', colnames (module_scores))
rownames (module_scores) <- gsub ('\\.', '-', rownames (module_scores) )
meta_data <- all_data@meta.data [match (colnames (module_scores), colnames (all_data) ), ]
module_signal <- CreateSeuratObject ( module_scores, meta.data = meta_data )
heat_signal <- ML$incorporate_TSC (module_signal)
select_signal <- heat_signal$select == 'select'

p5 <- DEG$seurat_heat_highlight (heat_signal, select_signal, rownames (heat_signal),
                                 c('new_cluster2'), average=T,
                                 heat_name=c('norm count', 'norm count'),
                                 col_legend_labels=c('cell type'),
                                 col_split = c(NA, 1),
                                 show_col_names = c(TRUE, FALSE),
                                 width = c(10, 10),normalize_data=T,
                                 col_rotation=c(90, 0), return_sep=T, 
                                 cluster_col=c(TRUE, TRUE))

p5_final <- ComplexHeatmap::draw (p5[[1]] + p5[[2]], 
                                  ht_gap = unit(c(0.3, 0.5), "cm"))

# ----------figure S3F: WGCNA modules----------
color_row <- read.csv ( paste ( sup_save_dir2, 'WGCNA/module_genes.csv' , sep='/'), row.names=1)
gene_list <- lapply (as.list (colnames (color_row) ), function (x) {unique (color_row [, x]) })
names (gene_list) <- WG$colors2labels (colnames (color_row), prefix='GC')

module_score <- KG$get_module_score (all_data, paste (sup_save_dir2, 
                'WGCNA/Data_module_score.csv', sep='/'), pgenes=gene_list)
colnames (module_score) <- gsub ('^X', '', colnames (module_score))
module_data <- CreateSeuratObject ( module_score, meta.data = all_data@meta.data )
heat_data <- ML$incorporate_TSC (module_data)
select_cells <- heat_data$select == 'select'

sel_genes <- rownames (module_data)
names (sel_genes) <- c('ECM', 'steroid', 'metabolism', 'ECM', 'NA', 'cytokine',
                       'immune', 'NA', 'immune', 'pluripotency',
                       'pluripotency')
gene_order <- c('metabolism', 'NA', 'pluripotency', 'ECM', 'immune', 'cytokine', 'steroid')

p6 <- DEG$seurat_heat_highlight (heat_data, select_cells, sel_genes,
                                 c('new_cluster2'), average=T,
                                 heat_name=c('norm count', 'norm count'),
                                 col_legend_labels=c('cell type'),
                                 row_legend_labels=c('WGCNA clusters'),
                                 width = c(10, 10),normalize_data=T,
                                 col_rotation=c(90, 0), return_sep=T, 
                                 row_reorder_levels = gene_order,
                                 col_split = c(NA, 1),
                                 show_col_names = c(TRUE, FALSE),
                                 cluster_col=c(TRUE, TRUE))
p6_final <- ComplexHeatmap::draw (p6[[1]] + p6[[2]], 
                                  ht_gap = unit(c(0.3, 0.5), "cm"))

# ----------figure S2G: cell cycle----------
CC <- modules::use ('cycle_analysis.R')
exp_mat <- as.matrix (GetAssayData (all_data, assay='RNA', slot='data'))
ans <- CC$get_phase_score (exp_mat [, all_data$new_cluster2 != 'hTSC'])
p7 <- CC$make_cycle_heat (ans, all_data, feature='new_cluster2')

# ----------merge everything----------
grob_list <- list (p1[[1]] + labs (fill='cell type'), 
                   p2, p3, p4, p5_final, p6_final, p7)
lay_mat <- matrix(c(1, 1, 2, 2, 3, 3, 
                    4, 4, 4, 4, 4, 4,
                    5, 5, 5, 5, 7, 7,
                    6 ,6 ,6 ,6 ,7 ,7 
                    ),
                  nrow=6) %>% t()
FO$arrange_plots (grob_list, paste (sup_save_dir, 'final_figureS4.pdf', sep='/'), lay_mat, 
                  plot_width=3.5, plot_height=7)
