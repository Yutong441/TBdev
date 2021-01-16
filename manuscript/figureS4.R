# generate the figure S4 of the manuscript
# Identity mapping of Trophoblast stem cells

setwd ('..')
devtools::load_all ()
library (tidyverse)
library (org.Hs.eg.db)
library (GOSemSim)

root_dir <- '/mnt/c/Users/Yutong/Documents/bioinformatics/reproduction/'
root <- paste (root_dir, 'results/', sep='/')
merge_dir <- paste (root, 'XLYBPZ_Dylan_dir', sep='/')
save_dir <- paste (root, 'manuscript/figure3', sep='/')
sup_save_dir3 <- paste (root, 'manuscript/figureS3', sep='/')
sup_save_dir <- paste (root, 'manuscript/figureS4', sep='/')
x <- load (paste (merge_dir, 'final_merged_tb.Robj', sep='/') )
all_data <- get (x)
TSC_cells <- all_data$revised %in% c('hTSC_OKAE', 'hTSC_TURCO') 

# ----------figure 3SA: diffusion map ----------
data (CT)
all_data2 <- all_data [, ! (all_data$broad_type %in% c(CT$non_emb_lineage,
                                                    CT$pre_imp_lineage))]
all_data2 <- run_dim_red (all_data2, run_diff_map=T, var_scale=T,
                          normalize=F, find_var_features=T, run_umap=F)
# filter out cells that interfere with aesthetics
show_data <- all_data2 [, !(all_data2$assigned_cluster == 'uCTB' ) ]
highlight <- show_data$assigned_cluster %in% CT$in_vitro_cells

devtools::load_all ()
p1 <- plot_dim_red (show_data, group.by= c('broad_type'), DR='diff_map' , return_sep=F,
                    size_highlight=highlight, highlight_font=2, dims=c(1,2),
                    nudge_ratio=0.2, move_x=0, move_y=1, further_repel=F,
                    repel_force=10, reverse_x=T, length_ratio=0.1)
p1

# ----------figure S4B: over-representation for Okae----------
selected_cells <- c('eTB', 'iTB', 'aTB', 'vCTB1', 'hTSC_OKAE')
d <- godata(org.Hs.eg.db, ont="BP")
all_data$new_cluster2 <- as.character (all_data$assigned_cluster)
OKAE_data <- all_data [, all_data$new_cluster2 %in% selected_cells]
OKAE_data$new_cluster2 <- partial_relevel (OKAE_data$new_cluster2)

markers <- find_DE_genes (OKAE_data, sup_save_dir3, feature='new_cluster2', label='OKAE')
kk_okae <- compare_cluster_enrichment (markers, d, org.Hs.eg.db, enrich_area='KEGG')
set.seed(100)
p2<- display_cluster_enrichment (kk_okae, show_graph='emap', feature_vec=
                                     markers$group, show_num=20) + labs (fill='')
p2

# ----------figure S4C: for Turco----------
selected_cells <- c('vCTB1', 'vCTB2', 'vCTB3','hTSC_TURCO')
TURCO_data <- all_data [, all_data$new_cluster2 %in% selected_cells]
TURCO_data$new_cluster2 <- partial_relevel (TURCO_data$new_cluster2)
markers <- find_DE_genes (TURCO_data, sup_save_dir3, feature='new_cluster2', label='TURCO')
kk_turco <- compare_cluster_enrichment (markers, d, org.Hs.eg.db, enrich_area='KEGG')
set.seed(100)
p3 <- display_cluster_enrichment (kk_turco, show_graph='emap', feature_vec=
                                     markers$group, show_num=20) + labs (fill='')

# ----------figure S4D: probability----------
sup_save_dir2 <- paste (root, 'manuscript/figureS2', sep='/')
BRGP_prob <- read.csv (paste (sup_save_dir2, 'result/BRGP_likelihood.csv', sep='/') )
show_meta <- all_data@meta.data [!is.na (all_data$MGP_PT) & !all_data$broad_type %in% c('EPI', 'PE'),]
select_cells <- c('eICM', 'aICM', 'eTB', 'iTB', 'aTB', 'hTSC_OKAE', 'hTSC_TURCO', 'cleavage',
                  'eEVT', 'aEVT', 'eCTB', 'vCTB1', 'vCTB2', 'vCTB3',
                  'hESC', 'hESC_YAN')

devtools::load_all ()
p4 <- plot_prob_line (BRGP_prob, select_cells, CT$in_vitro_cells,
                meta=show_meta[show_meta$broad_type!='STB',], 
                vjust=1e-4, thickness=1e-4, normalize_data=T)
select_cells <- c('cleavage', 'eICM', 'aICM', 'eTB', 'iTB', 'aTB', 'hTSC_OKAE', 'hTSC_TURCO',
                  'hESC', 'hESC_YAN', 'eCTB', 'sCTB', 'aSTB1', 'aSTB2', 'aSTB3')
p5 <- plot_prob_line (BRGP_prob, select_cells, CT$in_vitro_cells,
                meta=show_meta [show_meta$broad_type!='EVT',], 
                vjust=1e-4, thickness=1e-4, sel_branch='branch2', normalize_data=T)

# ----------figure S4E: pathway module----------
sup_save_dir2 <- paste (root, 'manuscript/figureS2', sep='/')
module_scores <- get_module_score (all_data, save_path=paste (sup_save_dir2, 
                                        'Data_module_scores.csv', sep='/'))
colnames (module_scores) <- gsub ('^X', '', colnames (module_scores))
rownames (module_scores) <- gsub ('\\.', '-', rownames (module_scores) )
meta_data <- all_data@meta.data [match (colnames (module_scores), colnames (all_data) ), ]
module_signal <- Seurat::CreateSeuratObject ( module_scores, meta.data = meta_data )
heat_signal <- incorporate_TSC (module_signal)
select_signal <- heat_signal$select == 'select'

seurat_param <- list (
        heat_name=c('norm count'),
        column_legend_labels=c('cell type'),
        main_width = c(10),
        main_height= 12,
        column_split = c(NA, 1),
        column_rotation=c(90, 0),
        show_column_names = c(T, F),
        show_column_anna=c(T, F),
        cluster_column=c(T, T),
        column_title_side=c('bottom', 'top'),
        grid_height=5,heat_grid_height=8
)
p6 <- seurat_heat_highlight (heat_signal, select_signal, rownames (heat_signal),
                             c('new_cluster2'), average=T, return_sep=T,
                             seurat_heat_params=seurat_param, row_scale=T)

p6_final <- ComplexHeatmap::draw (p6[[1]] + p6[[2]], 
                                  ht_gap = unit(c(0.3, 0.5), "cm"))

# ----------figure S3F: WGCNA modules----------
color_row <- read.csv ( paste ( sup_save_dir2, 'WGCNA/module_genes.csv' , sep='/'), row.names=1)
gene_list <- lapply (as.list (colnames (color_row) ), function (x) {unique (color_row [, x]) })
names (gene_list) <- colors2labels (colnames (color_row), prefix='GC')

module_score <- get_module_score (all_data, paste (sup_save_dir2, 
                'WGCNA/Data_module_score.csv', sep='/'), pgenes=gene_list)
colnames (module_score) <- gsub ('^X', '', colnames (module_score))
module_data <- Seurat::CreateSeuratObject ( module_score, meta.data = all_data@meta.data )
heat_data <- incorporate_TSC (module_data)
select_cells <- heat_data$select == 'select'

sel_genes <- rownames (module_data)
names (sel_genes) <- c('ECM', 'steroid', 'metabolism', 'ECM', 'NA', 'cytokine',
                       'immune', 'NA', 'immune', 'pluripotency',
                       'pluripotency')
gene_order <- c('metabolism', 'NA', 'pluripotency', 'ECM', 'immune', 'cytokine', 'steroid')

seurat_param2 <- list (
        heat_name=c('norm count'),
        column_legend_labels=c('gene cluster'),
        main_width = c(10),
        main_height= 12,
        column_split = c(NA, 1),
        column_rotation=c(90, 0),
        show_column_names = c(T, F),
        show_column_anna=c(T, F),
        cluster_column=c(T, T),
        column_title_side=c('bottom', 'top'),
        grid_height=5, heat_grid_height=12
)
p7 <- seurat_heat_highlight (heat_data, select_cells, sel_genes,
                             c('new_cluster2'), average=T, return_sep=T,
                             seurat_heat_params=seurat_param, row_scale=T)

p7_final <- ComplexHeatmap::draw (p7[[1]] + p7[[2]], 
                                  ht_gap = unit(c(0.3, 0.5), "cm"))

# ----------figure S2G: cell cycle----------
exp_mat <- as.matrix (Seurat::GetAssayData (all_data, assay='RNA', slot='data'))
ans <- get_phase_score (exp_mat)
p8 <- make_cycle_heat (ans, all_data, feature='new_cluster2',automatic=F, grid_height=5)

# ----------merge everything----------
grob_list <- list (p1, p2, p3, 
                   p4 +theme (aspect.ratio=0.7) + labs (fill='cell type', color='cell type'), 
                   p5+theme (aspect.ratio=0.7)+ labs (fill='cell type', color='cell type'), 
                   p6_final, p7_final, p8)
lay_mat <- matrix(c(1, 1, 2, 2, 3, 3, 
                    4, 4, 4, 5, 5, 5,
                    6, 6, 6, 6, 8, 8,
                    7 ,7 ,7 ,7 ,8 ,8 
                    ),
                  nrow=6) %>% t()
arrange_plots (grob_list, paste (sup_save_dir, 'final_figureS4.pdf', sep='/'), lay_mat, 
                  plot_width=3.5, plot_height=7)
