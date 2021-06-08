devtools::load_all ('..', export_all=T)
library (tidyverse)
root_dir <- '/mnt/c/Users/Yutong/Documents/bioinformatics/reproduction/'
root <- paste (root_dir, 'results/', sep='/')
merge_dir <- paste (root, 'XLYBPZ_Dylan_dir', sep='/')
save_dir <- paste (root, 'manuscript/ISSCR_poster', sep='/')
save_dir1 <- paste (root, 'manuscript/figure1', sep='/')
save_dir4 <- paste (root, 'manuscript/figure4', sep='/')

all_data <- readRDS(paste (merge_dir, 'merged_blastoid.rds', sep='/') )
plotlist <- list ()

# ----------change labels----------
# exclude irrelevant cells
all_data <- all_data [, !all_data$broad_type %in% c('hTSC-TURCO', 'bPE-YANA')]
all_data$broad_type <- gsub ('hESC', 'hPSC', all_data$broad_type)
all_data$final_cluster <- gsub ('hESC', 'hPSC', all_data$final_cluster)
all_data$broad_type <- gsub ('^TB$', 'TE', all_data$broad_type)
all_data$final_cluster <- gsub ('^TB$', 'TE', all_data$final_cluster)

data (format_conf)
ori_order <- gsub ('hESC', 'hPSC', format_conf$cell_order)
ori_order <- gsub ('^TB$', 'TE', ori_order)
new_order <- c(ori_order, 'hPSC', 'nPSC', 'bEPI-YANA', 'hTSC-OKAE','bTSC',
               'bTB-YANA')

all_data$broad_type <- TBdev::partial_relevel (all_data$broad_type, new_order)
all_data$final_cluster <- TBdev::partial_relevel (all_data$final_cluster, new_order)
TB_data <- all_data [, all_data$date != 'in_vitro']

# ----------Aesthetic settings----------
blastoid_color <- c('nESC'='#71d300', 'bTB-YANA'='#ff1900', 'bTSC'='#ff8500', 'bEPI-YANA'='#00ff9d')
AP <- list (color_vec=c(format_conf$color_vec, blastoid_color),
            cell_order=new_order, edge_stroke=1e-5, 
            point_edge_color=alpha ('gray', 0.), fontsize=22,
            point_fontsize=22/3)
AP$color_vec <- AP$color_vec [order (partial_relevel (names ( AP$color_vec),
                                                      AP$cell_order))]

# ----------Title logo----------
plotlist [['cam_logo']] <- jpeg::readJPEG (paste (save_dir, 'cambridge_logo.jpg', sep='/'))
plotlist [['pdn_logo']] <- png::readPNG (paste (save_dir, 'PDN_logo_white.png', sep='/'))

# This script generates poster for the PDN symposium
# ----------Poster figure 1: Introduction----------
plotlist [['introduction_cell']] <- png::readPNG (paste (save_dir, 
                                'intro_fig.png', sep='/'))

# ----------Poster figure 2: integration----------
plot_dim_red (TB_data, group.by= c('broad_type', 'date'), DR='pca' , dims=c(1,2),
              return_sep=T, nudge_ratio=0., plot_type='dim_red_sim',
              seg_color='black', nudge_ortho=0.7, 
              length_ratio=0.1, AP=AP, further_repel=F) -> p1

plotlist [['integrate_type']] <- p1[[1]] + labs (fill='')
plotlist [['integrate_date']] <- p1[[2]] + labs (fill='')

data (lineage_markers)
show_genes <- lineage_markers [!names (lineage_markers) %in% c('STR', 'PE')]
show_genes <- show_genes [!show_genes %in% c('PRDM14', 'ELF5', 'ETS2', 'HAND1',
                                             'MMP9', 'KLF17', 'CGB2', 'LHB')]
seurat_heat (TB_data, color_row=show_genes, group.by = c('broad_type'), 
             slot_data='data', heat_name='norm count', center_scale=T,
             column_legend_labels=c('cell type'), row_scale=T, 
             main_width=24, main_height=16, automatic=F, 
             row_legend_labels ='lineage gene',
             column_rotation=90, AP=AP) -> plotlist [['integrate_heat']]

# ----------Poster figure 4: Trajectory inference----------
# GPLVM + trajectory
sup_save_dir2 <- paste (root, 'manuscript/figureS2', sep='/')
epg <- read.csv (paste (sup_save_dir2, 'result/STREAM_graph.csv', sep='/'))
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
dim_red_3D_traj (plot_met, 'PT1', 'PT2', 'PT3', 'broad_type', epg, 'x', 'y',
                    'z', 'branch_name', all_theta=50, all_phi=0, further_repel=T,
                    repel_force=0.5, lab_just=c(0.08, 0.02, 0.02), magnify_text=1.3, 
                    label_traj_text=label_epg, hor_just=0.1, dim_vjust=4, AP=AP) + 
labs (fill='cell type') +theme (legend.position='none', aspect.ratio=0.6) -> plotlist [['traj_gp']]

# signaling pathway
module_scores <- get_module_score (all_data,
                              save_path=paste (sup_save_dir2, 'Data_module_scores.csv', sep='/'))
rownames (metadata) <- gsub ('\\.[0-9]+\\.[0-9]+$', '', rownames (metadata) )
module_scores <- module_scores [, match (rownames (metadata), colnames (module_scores) ) ]
rownames (module_scores) <- gsub ('\\.', '-', rownames (module_scores) )
path_seurat <- Seurat::CreateSeuratObject (module_scores, meta.data=metadata)

data (format_conf)
exclude_path <- c('TCA', 'apoptosis', 'TJ', 'adherens', 'gap', 'JAK-STAT')
sel_path <- rownames (path_seurat) [!rownames (path_seurat) %in% exclude_path]

seurat_heat (path_seurat, group.by=c('epil_branch','broad_type'),
             color_row=sel_path, row_scale=T,
             column_reorder_levels = list (format_conf$branch_order, format_conf$cell_order),
             column_legend_labels= c('branch', 'cell type'), 
             heat_name = 'module score', center_scale=T,
             group_order = order (path_seurat$MGP_PT), automatic=F,
             show_column_bars = c(T,T), cluster_rows=T, AP=AP
) -> plotlist [['traj_signal']]

# ----------figure 5: compare in vitro and in vivo ----------
# gene expression pattern
in_vitro_cells <- c('hESC', 'hESC-YAN', 'nESC', 'bEPI-YANA', 'hTSC-OKAE', 'bTSC', 'bTB-YANA')
in_vitro <- all_data [, all_data$broad_type %in% in_vitro_cells ]
plot_genes <- c('POU5F1', 'SOX2', 'NANOG', 'KRT7', 'KRT18', 'TFAP2C',
                'TFAP2A','GATA3','GATA2', 'ELF5', 'CDX2', 'TEAD4', 'HLA-A',
                'HLA-B')
seurat_violin (in_vitro, features=plot_genes, group.by='broad_type',
               num_col=7, plot_type='box', free_xy='free_y', lower_b=0, AP=AP)+
             theme (aspect.ratio=1.2, legend.position='none')-> plotlist [['hTSC_gene']]

# load fitting probability
fit_chance <- read.csv (paste (root_dir, 'results/Yu_2021/fit_chance.csv', sep='/'))
fit_chance2 <- read.csv (paste (root_dir, 'results/Yanagida_2021/fit_chance.csv', sep='/'))
fit_chance2 %>% select (!X) %>% cbind (fit_chance) %>% 
        column_to_rownames (var='X') -> fit_chance
show_meta <- cbind (metadata [, !colnames (metadata) %in% colnames (fit_chance)], 
                    fit_chance [match (rownames (metadata), rownames (fit_chance)),])

# in GPLVM
show_meta %>% slice_min (PT3, n=nrow(show_meta)-3) %>% 
        drop_na (hTSC_OKAE) %>% 
        normalize_prob(c('hESC', 'hTSC.OKAE', 'bTB.YANA', 'bTSC')) %>%
        rename('hTSC-OKAE'='hTSC.OKAE') %>% rename ('bTB-YANA'='bTB.YANA') %>%
        dim_red_3D_traj ('PT1', 'PT2', 'PT3', c('hESC', 'hTSC-OKAE', 'bTSC', 'bTB-YANA'), 
                         epg, 'x', 'y', 'z', 'branch_name', all_theta=50, all_phi=0, 
                         show_label=T, further_repel=F, 
                         lab_just=c(0.08, 0.02, 0.02), label_col='broad_type',
                         label_traj_text=label_epg, num_col=2, AP=AP,
                         hor_just=0.1, dim_vjust=3)+ coord_cartesian (clip='off')+
        labs (fill='relative probability') +
        scale_color_manual (values=rep('black',3))+guides(color=F) -> plotlist [['hTSC_gp']]

# ----------TF network----------
# GRN inference
color_row <- read.csv ( paste ( save_dir4,  'WGCNA/module_genes.csv' , sep='/'), row.names=1)
gene_list <- lapply (as.list (colnames (color_row) ), function (x) {
                             unique (color_row [, x]) })
names (gene_list) <- colnames (color_row)
plot_gene <- list (TB=gene_list$GC1, CTB=gene_list$GC6, STB=gene_list$GC3, EVT=gene_list$GC10)
markers <- find_DE_genes (all_data, save_dir1, group.by='broad_type', label='all_vivo')

data (TF)
tf_data <- all_data [TF, all_data$date!='in_vitro']
pnet <- custom_net_diff_nets (tf_data, plot_gene, markers, nudge_ratio=0.3, 
                              size_thres=0.2, AP=AP, ranges=c(0, 30))
pnet <- lapply (pnet, function (gx){gx+coord_cartesian (clip='off')})
plotlist [['traj_network']] <- ggpubr::ggarrange (plotlist=pnet[2:4], nrow=1, 
                                                  common.legend=T, legend='right')

# ----------siRNA screen----------
img_dir <- paste (root_dir, 'results/imaging/clono', sep='/')
proc_dat <- read.csv (paste (img_dir, 'combined_proc.csv', sep='/'))
proc_dat %>% group_by (condition, Date) %>%
        summarise (total_num = sum(rel_num_colony)) -> by_run

cell_lab <- read.csv('../inst/python/clono/TF_cell.csv')
by_run$celltype <- cell_lab$celltype [match (by_run$condition, toupper (cell_lab$TF))]
# filter out poor quality replicates:
by_run %>% filter (! Date %in% c('run3') ) -> by_run

cluster_level <- c('GFP', 'PITX1', 'CTNNB1', 'TFAP2C', 'HMOX1', 'ZNF750',
                   'SP6', 'MSX2', 'OVOL1', 'GRHL1', 'GCM1', 'TEAD3', 'GATA2',
                   'ZFHX3', 'GATA3', 'IRX4', 'NFE2L3', 'NR2F2', 'TFEB',
                   'TRIM28', 'MAZ', 'SSRP1', 'PCBP2', 'HMGA1')
duals <- c('GATA2_CEBPA', 'GATA2_GATA3', 'MSX2_NR2F2', 'MSX2_OVOL1', 'TEAD3_GCM1')
by_run %>% filter (!Date %in% c('run6', 'run8') & !condition %in% duals ) %>%
plot_colony_count (y_val='total_num', plot_type='box', paired_test=T, 
                   new_level=cluster_level,
                   plot_pval='p.signif', display_shape=F)+
ylab ('normalized number of colonies')-> plotlist [['single_siRNA']]

by_run %>% filter (condition %in% c(duals, 'GATA2', 'GATA3', 'MSX2', 'OVOL1',
                                    'TEAD3', 'GCM1', 'GFP') ) %>%
plot_colony_count (y_val='total_num', plot_type='box', paired_test=T,
                   plot_pval='p.signif', display_shape=F, new_level=cluster_level)+
theme (aspect.ratio=1)-> plotlist [['dual_siRNA']]

# graphical abstract
plotlist [['abstract']] <- png::readPNG (paste (save_dir, 'graphical_abstract.png', sep='/'))

# ----------combine all plots----------
poster_text <- rjson::fromJSON(file = paste (save_dir, 'poster_text.json', sep='/') )
plotlist_t <- c(plotlist, poster_text)
config_arr <- read.csv (paste (save_dir, 'poster_config.csv', sep='/'))
arrange_poster (plotlist_t, paste (save_dir, 'ISSCR_poster.pdf', sep='/'), config_arr, 
                plot_width=7.5/6, plot_height=7/24)
