devtools::load_all('..', export_all=F)
library (tidyverse)
library (Seurat)

root_dir <- '/mnt/c/Users/Yutong/Documents/bioinformatics/reproduction/'
root <- paste (root_dir, 'results/', sep='/')
merge_dir <- paste (root, 'XLYBPZ_Dylan_dir', sep='/')
save_dir <- paste (root, 'manuscript/figure3', sep='/')
sup_save_dir <- paste (root, 'manuscript/figureS3', sep='/')
all_data <- readRDS (file=paste (merge_dir, 'merged_blastoid.rds', sep='/'))

# ----------aesthetics----------
data (format_conf)
new_order <- c(format_conf$cell_order, 'hESC', 'nESC', 'bEPI-YANA', 'hTSC-OKAE','bTSC', 'bTB-YANA')
blastoid_color <- c('nESC'='#71d300', 'bTB-YANA'='#ff1900', 'bTSC'='#ff8500', 'bEPI-YANA'='#00ff9d')
new_color <- list (color_vec=c(format_conf$color_vec, blastoid_color), cell_order=new_order)
new_color$color_vec <- new_color$color_vec [order (partial_relevel (names (
                                        new_color$color_vec), new_color$cell_order))]

# ----------figure 3A ----------
in_vitro_cells <- c('hESC', 'hESC-YAN', 'nESC', 'bEPI-YANA', 'hTSC-OKAE', 'bTSC', 'bTB-YANA')
in_vitro <- all_data [, all_data$broad_type %in% in_vitro_cells ]
plot_genes <- c('POU5F1', 'SOX2', 'NANOG', 
                'TDGF1', 'PRDM14', 'GDF3', 'KRT7', 'GATA3','GATA2', 'TFAP2C', 'CDX2', 'TEAD4', 
                'ENPEP', 'TACSTD2','SIGLEC6', 'ISL1', 'GABRP', 'VTCN1',
                'HLA-A', 'HLA-B', 'HLA-C'
)

p1 <- seurat_violin (in_vitro, features=plot_genes, group.by='final_cluster',
                     num_col=3, free_xy='free_y', lower_b=0, plot_type='box', AP=new_color)

# ----------figure 3B ----------
data (CT)
# do not show the blastoid cells which will obscure Turco and Okae cells
pca_data <- all_data [, !all_data$broad_type %in% c('nESC', 'bPE-YANA', 'hTSC-TURCO')]
dim_red_dat <- pca_data [, !pca_data$broad_type %in% c('bTSC', 'bTB-YANA')]
p2 <- plot_dim_red (dim_red_dat, group.by= c('broad_type'), DR='pca', return_sep=T, 
                    size_highlight=dim_red_dat$broad_type=='hTSC-OKAE', 
                    highlight_ratio=2, dims=c(1,2), nudge_ratio=0.15,
                    plot_type='dim_red_sim', seg_color='black', 
                    AP=c(list (point_edge_color='gray'), new_color))

dim_red_dat <- pca_data [, !pca_data$broad_type %in% c('hTSC-OKAE', 'bTB-YANA')]
p2_2 <- plot_dim_red (dim_red_dat, group.by= c('broad_type'), DR='pca', return_sep=T, 
                    size_highlight=dim_red_dat$broad_type=='bTSC', 
                    highlight_ratio=2, dims=c(1,2), nudge_ratio=0.15,
                    plot_type='dim_red_sim', seg_color='black', 
                    AP=c(list (point_edge_color='gray'), new_color))

dim_red_dat <- pca_data [, !pca_data$broad_type %in% c('hTSC-OKAE', 'bTSC')]
p2_3 <- plot_dim_red (dim_red_dat, group.by= c('broad_type'), DR='pca', return_sep=T, 
                    size_highlight=dim_red_dat$broad_type=='bTB-YANA', 
                    highlight_ratio=2, dims=c(1,2), nudge_ratio=0.15,
                    plot_type='dim_red_sim', seg_color='black', 
                    AP=c(list (point_edge_color='gray'), new_color))

# ----------figure 3C----------
sup_save_dir2 <- paste (root, 'manuscript/figureS2', sep='/')
epg <- read.csv (paste (sup_save_dir2, 'result/STREAM_graph.csv', sep='/'))
show_meta <- all_data@meta.data
show_meta %>% filter (!is.na (MGP_PT) ) %>% filter (!broad_type %in% c('PE', 'EPI') ) -> show_meta

new_name <- c('TB_stem', 'EVT_branch', 'STB_branch', NA)
old_name <- c('S2,S1', 'S0,S1', 'S3,S1', 'blank')
data (format_conf)
epg$branch_name <- partial_relevel (new_name[match (epg$branch, old_name )], format_conf$branch_order)

epg %>% arrange (branch_name) %>% mutate( index = 1:nrow (epg)) %>% 
        group_by (branch) %>% slice_min (index, n=1) %>% drop_na() %>% 
        as.data.frame () %>% dplyr::select (x, y, z, branch_name) %>% 
        magrittr::set_colnames (c('x', 'y', 'z', 'feature')) -> label_epg

fit_chance <- read.csv (paste (root_dir, 'results/Yu_2021/fit_chance.csv', sep='/'))
fit_chance2 <- read.csv (paste (root_dir, 'results/Yanagida_2021/fit_chance.csv', sep='/'))

fit_chance2 %>% select (!X) %>% cbind (fit_chance) %>% 
        column_to_rownames (var='X') -> fit_chance

rownames (show_meta) <- gsub ('\\.[0-9]+$', '', rownames (show_meta))
show_meta <- cbind (show_meta [, !colnames (show_meta) %in% colnames (fit_chance)], 
                    fit_chance [match (rownames (show_meta), rownames (fit_chance)),])

show_meta %>% slice_min (PT3, n=nrow(show_meta)-3) %>% 
        normalize_prob(c('hTSC.OKAE', 'bTB.YANA', 'bTSC')) %>%
        rename('hTSC-OKAE'='hTSC.OKAE') %>% rename ('bTB-YANA'='bTB.YANA') %>%
        dim_red_3D_traj ('PT1', 'PT2', 'PT3', c('hTSC-OKAE', 'bTSC', 'bTB-YANA'), epg, 'x', 'y',
        'z', 'branch_name', all_theta=50, all_phi=0, show_label=T, further_repel=T,
        repel_force=0.2, lab_just=c(0.08, 0.02, 0.02), label_col='broad_type',
        label_traj_text=label_epg, num_col=1, AP=c(list (point_edge_color='white'), new_color), 
        hor_just=0.1, dim_vjust=3) + 
        labs (fill='relative probability') +
        scale_color_manual (values=rep('black',3))+guides(color=F) -> p3

# ----------figure 3D: cell-cell correlation ----------
in_vitro_cells <- c('hESC', 'nESC', 'bEPI-YANA', 'hTSC-OKAE', 'bTSC', 'bTB-YANA')
select_cells2 <- all_data$broad_type %in% in_vitro_cells
select_cells1 <- all_data$date != 'in_vitro'
#all_cor <- compute_all_cor (all_data, method= 'correlation', assay='RNA',
#                            select_cells2=select_cells2,
#                            select_cells1=select_cells1)
#saveRDS(all_cor, paste (root_dir, 'results/Yu_2021/cell_cor.rds', sep='/'))
all_cor <- readRDS (paste (root_dir, 'results/Yu_2021/cell_cor.rds', sep='/'))
p4 <- cell_violin (all_cor, all_data@meta.data, c('final_cluster', 'final_cluster'), 
                   box_plot=T, num_col=2, legend_col=1, column_scale=T, AP=new_color)

# ----------arrange all figures----------
grob_list <- list (p1 + theme (aspect.ratio=1., 
                        axis.title.x=element_blank(), legend.position='top') + labs (fill =''), 
                   p2[[1]] +labs (fill=''), 
                   p2_2[[1]] +labs (fill=''), 
                   p2_3[[1]] +labs (fill=''), 
                   ggplot() + theme_minimal (),
                   p3+ theme(legend.position='top'), 
                   p4+ labs (fill='') + theme (axis.title.x=element_blank(), aspect.ratio=0.8)
)
lay_mat <- matrix(c(
                    1, 1, 1, 2, 2, 2, 
                    1, 1, 1, 3, 3, 3, 
                    4, 4, 4, 5, 5, 5,
                    6, 6, 7, 7, 7, 7,
                    6 ,6 ,7 ,7 ,7 ,7
                    ),
                  nrow=6) %>% t()
arrange_plots (grob_list, paste (save_dir, 'final_figure3.pdf', sep='/'),
               lay_mat, plot_width=3.2, plot_height=8)
save_indiv_plots (grob_list, paste (save_dir, 'figure3', sep='/'), lay_mat,
                  plot_width=3.2, plot_height=8)
