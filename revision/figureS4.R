# generate the figure S4 of the manuscript
# Identity mapping of Trophoblast stem cells

setwd ('..')
#roxygen2::roxygenise()
devtools::load_all ('.', export_all=F)
library (tidyverse)
library (org.Hs.eg.db)
library (Seurat)

root_dir <- '/mnt/c/Users/Yutong/Documents/bioinformatics/reproduction/'
root <- paste (root_dir, 'results/', sep='/')
merge_dir <- paste (root, 'XLYBPZ_Dylan_dir', sep='/')
save_dir <- paste (root, 'manuscript/figure3', sep='/')
sup_save_dir <- paste (root, 'manuscript/figureS4', sep='/')
all_data_sub <- get(load(paste (merge_dir, 'final_merged_tb.Robj', sep='/') ))
all_data_sub <- all_data_sub [, all_data_sub$final_cluster!='uCTB']
all_data <- readRDS(paste (merge_dir, 'merged_blastoid.rds', sep='/') )

# ----------aesthetics----------
data (format_conf)
new_order <- c(format_conf$cell_order, 'hESC', 'nESC', 'bEPI-YANA', 'hTSC-OKAE','bTSC', 'bTB-YANA')
blastoid_color <- c('nESC'='#71d300', 'bTB-YANA'='#ff1900', 'bTSC'='#ff8500', 'bEPI-YANA'='#00ff9d')
new_color <- list (color_vec=c(format_conf$color_vec, blastoid_color), cell_order=new_order)
new_color$color_vec <- new_color$color_vec [order (partial_relevel (names (
                                        new_color$color_vec), new_color$cell_order))]

# ----------figure 3SA: diffusion map ----------
data (CT)
no_btsc <- all_data [, !all_data$broad_type %in% c('nESC', 'hTSC-TURCO', 'bPE-YANA', 
                                                   'hESC', 'hESC-YAN', 'bEPI-YANA')]
all_data2 <- run_dim_red (no_btsc, run_diff_map=T, var_scale=T,
                          normalize=F, find_var_features=T, run_umap=F, 
                          select_cells=!no_btsc$date %in% 'in_vitro')
# remove irrelevant cells
show_data <- all_data2 [, !(all_data2$broad_type %in% c(
                         CT$non_emb_lineage, CT$pre_imp_lineage))]

dim_red_dat <- show_data [, !show_data$broad_type %in% c('bTSC', 'bTB-YANA')]
p1 <- plot_dim_red (dim_red_dat, group.by= c('broad_type'), DR='diff_map' , return_sep=T,
                    size_highlight=dim_red_dat$broad_type=='hTSC-OKAE', 
                    dims=c(1,2), nudge_ratio=0.1, move_x=0, move_y=0, highlight_ratio=2,
                    further_repel=F, repel_force=10, reverse_x=T,
                    length_ratio=0.05, plot_type='dim_red_sim', AP=new_color)

dim_red_dat <- show_data [, !show_data$broad_type %in% c('hTSC-OKAE', 'bTB-YANA')]
p1_2 <- plot_dim_red (dim_red_dat, group.by= c('broad_type'), DR='diff_map' , return_sep=T,
                    size_highlight=dim_red_dat$broad_type=='bTSC', highlight_ratio=2,
                    dims=c(1,2), nudge_ratio=0.1, move_x=0, move_y=0,
                    further_repel=F, repel_force=10, reverse_x=T,
                    length_ratio=0.05, plot_type='dim_red_sim', AP=new_color)

dim_red_dat <- show_data [, !show_data$broad_type %in% c('bTSC', 'hTSC-OKAE')]
p1_3 <- plot_dim_red (dim_red_dat, group.by= c('broad_type'), DR='diff_map' , return_sep=T,
                    size_highlight=dim_red_dat$broad_type=='bTB-YANA', highlight_ratio=2,
                    dims=c(1,2), nudge_ratio=0.1, move_x=0, move_y=0,
                    further_repel=F, repel_force=10, reverse_x=T,
                    length_ratio=0.05, plot_type='dim_red_sim', AP=new_color)
rm (all_data2)

# ----------figure S4B-C: probability----------
in_vitro <- c( 'hTSC_OKAE', 'hTSC_TURCO', 'hESC', 'hESC_YAN')
in_vivo <- c('ICM', 'TB', 'CTB', 'STB', 'EVT')
sup_save_dir2 <- paste (root, 'manuscript/figureS2', sep='/')
BRGP_prob <- read.csv (paste (sup_save_dir2, 'result/cell_likelihood_broad_type.csv', sep='/'), row.names=1)
BRGP_prob %>% dplyr::select (all_of (c(in_vivo, 'hESC'))) %>%
        mutate_if (is.numeric, function(x){x/max(x)}) %>%
        magrittr::set_rownames (rownames (BRGP_prob)) -> BRGP_prob

sup_save_dir2 <- paste (root, 'manuscript/figureS2', sep='/')
epg <- read.csv (paste (sup_save_dir2, 'result/STREAM_graph.csv', sep='/'))
new_name <- c('TB_stem', 'EVT_branch', 'STB_branch', NA)
old_name <- c('S2,S1', 'S0,S1', 'S3,S1', 'blank')
data (format_conf)
epg$branch_name <- partial_relevel (new_name[match (epg$branch, old_name )], format_conf$branch_order)
epg %>% arrange (branch_name) %>% mutate( index = 1:nrow (epg)) %>% 
        group_by (branch) %>% slice_min (index, n=1) %>% drop_na() %>% 
        as.data.frame () %>% dplyr::select (x, y, z, branch_name) %>% 
        magrittr::set_colnames (c('x', 'y', 'z', 'feature')) -> label_epg

all_data_sub@meta.data [match (rownames (BRGP_prob), colnames (all_data_sub)),] %>%
        dplyr::select (!one_of (colnames (BRGP_prob))) %>% cbind (BRGP_prob) %>%
        filter (!broad_type %in% c('PE', 'EPI') ) -> show_meta

show_meta %>% slice_min (PT3, n=nrow(show_meta)-3) -> plot_met
p2 <- dim_red_3D_traj (plot_met, 'PT1', 'PT2', 'PT3', c(in_vivo, 'hESC'), epg, 'x', 'y',
                    'z', 'branch_name', all_theta=50, all_phi=0, further_repel=T,
                    repel_force=0.5, lab_just=c(0.08, 0.02, 0.02), hor_just=0.1, magnify_text=1.3, 
                    label_traj_text=label_epg, AP=list (point_edge_color='white'), num_col=2, dim_vjust=1.8) +
                    scale_color_manual (values=rep('black',3))+guides(color=F)+
                    labs (fill='relative probability')

# ----------figure S4C: specific subtypes----------
invivo <- all_data_sub [, all_data_sub$date!='in_vitro']
p3 <- plot_dim_red (invivo, group.by= c('final_cluster'), DR='pca', return_sep=T,
                    nudge_ratio=0.2, plot_type='dim_red_sim', seg_color='black')

# ----------figure S4D: invitro vs TB----------
save_dir3 <- paste (root, 'manuscript/figure3', sep='/')
invitro_cell <- c('hTSC-OKAE', 'bTSC', 'bTB-YANA')
invivo_comp <- c('TB', 'CTB', 'STB', 'EVT')
sel_data <- all_data [, all_data$broad_type %in% c(invitro_cell, invivo_comp)]
markers <- find_DE_genes (sel_data, save_dir3, group.by='broad_type', 
                          label='pairwise_blastoid', method='pairwise')

gsea <- run_GSEA_pairwise (markers, org.Hs.eg.db, paste (sup_save_dir, 'PSEA_pairwise.csv', sep='/'),
                   group1=invitro_cell, group2=invivo_comp)

devtools::load_all ('.', export_all=F)
rich_forest (gsea, markers, org.Hs.eg.db, show_num=2, show_gene_labels=4,
             shrink_ratio=0.7, band_width=0.3, band_ratio=3, extend_x_neg=0.4, AP=new_color)+
theme (aspect.ratio=1.4)-> p4
p4

# ----------merge everything----------
grob_list <- list (p1[[1]]+labs (fill=''), 
                   p1_2[[1]]+labs (fill=''),
                   p1_3[[1]]+labs (fill=''),
                   p2+theme (legend.position='top'), 
                   p3[[1]]+labs (fill=''), p4
)
lay_mat <- matrix(c(1,4, 
                    1,4,
                    1,4,
                    2,4,
                    2,4,
                    2,4,
                    3,5,
                    3,5,
                    3,5,
                    6,6,
                    6,6,
                    6,6,
                    6,6
                    ),
                  nrow=2) %>% t()
arrange_plots (grob_list, paste (sup_save_dir, 'final_figureS3.pdf', sep='/'), lay_mat, 
                  plot_width=9, plot_height=2.54, margin_width=1)
save_indiv_plots (grob_list, paste (sup_save_dir, 'figureS3', sep='/'),
                  lay_mat, plot_width=9, plot_height=2.54)
