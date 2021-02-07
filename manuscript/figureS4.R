# generate the figure S4 of the manuscript
# Identity mapping of Trophoblast stem cells

devtools::load_all ('..', export_all=F)
library (tidyverse)
library (org.Hs.eg.db)
library (GOSemSim)

root_dir <- '/mnt/c/Users/Yutong/Documents/bioinformatics/reproduction/'
root <- paste (root_dir, 'results/', sep='/')
merge_dir <- paste (root, 'XLYBPZ_Dylan_dir', sep='/')
save_dir <- paste (root, 'manuscript/figure3', sep='/')
sup_save_dir <- paste (root, 'manuscript/figureS4', sep='/')
all_data <- get (load (paste (merge_dir, 'final_merged_tb.Robj', sep='/') ))
all_data_sub <- all_data [, all_data$assigned_cluster!='uCTB']
invivo <- all_data_sub [, all_data_sub$date!='in_vitro']

# ----------figure 3SA: diffusion map ----------
data (CT)
all_data2 <- all_data [, ! (all_data$broad_type %in% c(
                        CT$non_emb_lineage, CT$pre_imp_lineage))]
all_data2 <- run_dim_red (all_data2, run_diff_map=T, var_scale=T,
                          normalize=F, find_var_features=T, run_umap=F)
show_data <- all_data2 [, !(all_data2$assigned_cluster == 'uCTB' ) ]
highlight <- show_data$final_cluster %in% CT$in_vitro_cells

p1 <- plot_dim_red (show_data, group.by= c('broad_type'), DR='diff_map' , return_sep=T,
                    size_highlight=highlight, dims=1:2,
                    nudge_ratio=0.2, move_x=0, move_y=0, further_repel=F,
                    repel_force=10, reverse_x=T, length_ratio=0.05, 
                    plot_type='dim_red_sim')
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

all_data@meta.data [match (rownames (BRGP_prob), colnames (all_data)),] %>%
        dplyr::select (!one_of (colnames (BRGP_prob))) %>% cbind (BRGP_prob) %>%
        filter (!broad_type %in% c('PE', 'EPI') ) -> show_meta

show_meta %>% slice_min (PT3, n=nrow(show_meta)-3) -> plot_met
devtools::load_all ('..', export_all=F)
p2 <- dim_red_3D_traj (plot_met, 'PT1', 'PT2', 'PT3', c(in_vivo, 'hESC'), epg, 'x', 'y',
                    'z', 'branch_name', all_theta=50, all_phi=0, further_repel=T,
                    repel_force=0.5, lab_just=c(0.08, 0.02, 0.02), hor_just=0.1, magnify_text=1.3, 
                    label_traj_text=label_epg, AP=list (point_edge_color='white'), num_col=2, dim_vjust=1.8) +
                    scale_color_manual (values=rep('black',3))+guides(color=F)+
                    labs (fill='relative probability')

# ----------figure S4D: specific subtypes----------
p3 <- plot_dim_red (invivo, group.by= c('final_cluster'), DR='pca', return_sep=T,
                    nudge_ratio=0.2, plot_type='dim_red_sim', seg_color='black')

# ----------figure S2G: cell cycle----------
exp_mat <- as.matrix (Seurat::GetAssayData (all_data_sub, assay='RNA', slot='data'))
ans <- get_phase_score (exp_mat)
p4 <- make_cycle_heat (ans, all_data, group.by='final_cluster',automatic=T,
                       grid_height=5, heat_grid_height=9,
                       direction='horizontal', main_height=9,
                       column_title_side='bottom')

# ----------merge everything----------
grob_list <- list (p1[[1]]+labs (fill=''), 
                   p2+theme (legend.position='top'), 
                   p3[[1]]+labs (fill=''), p4
)
lay_mat <- matrix(c(1, 2, 
                    3, 2,
                    4 ,4
                    ),
                  nrow=2) %>% t()
arrange_plots (grob_list, paste (sup_save_dir, 'final_figureS3.pdf', sep='/'), lay_mat, 
                  plot_width=9, plot_height=7)
save_indiv_plots (grob_list, paste (sup_save_dir, 'figureS3', sep='/'),
                  lay_mat, plot_width=9, plot_height=7)

# table S4
data (CT)
save_dir1 <- paste (root, 'manuscript/figure1', sep='/')
show_meta %>% dplyr::select ( all_of (c(in_vivo, in_vitro, 'final_cluster') ) ) %>%
        gather ( 'cell_type', 'prob', -final_cluster) %>%
        group_by (final_cluster, cell_type ) %>%
        summarise (max_prob = max (prob) ) %>%
        magrittr::set_colnames (c('reference', 'compare', 'max_prob')) %>%
        as.data.frame () %>%
        mutate (compare = partial_relevel (compare, CT$cell_order) ) %>%
        spread (compare, max_prob) %>%
        write.csv (paste (save_dir1, 'TableS4.csv', sep='/'), row.names=F)
