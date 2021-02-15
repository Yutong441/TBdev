devtools::load_all ('..', export_all=F)
library (tidyverse)
root_dir <- '/mnt/c/Users/Yutong/Documents/bioinformatics/reproduction/'
root <- paste (root_dir, 'results/', sep='/')
merge_dir <- paste (root, 'XLYBPZ_Dylan_dir', sep='/')
save_dir <- paste (root, 'manuscript/poster', sep='/')
save_dir1 <- paste (root, 'manuscript/figure1', sep='/')
save_dir4 <- paste (root, 'manuscript/figure4', sep='/')

all_data <- get(load (paste (merge_dir, 'final_merged_tb.Robj', sep='/') ))
TB_data <- all_data [, !c(all_data$revised %in% CT$in_vitro_cells)]
plotlist <- list ()

# ----------Title logo----------
plotlist [['cam_logo']] <- jpeg::readJPEG (paste (save_dir, 'cambridge_logo.jpg', sep='/'))
plotlist [['pdn_logo']] <- png::readPNG (paste (save_dir, 'PDN_logo_white.png', sep='/'))

# This script generates poster for the PDN symposium
# ----------Poster figure 1: Introduction----------
plotlist [['introduction_cell']] <- png::readPNG (paste (save_dir, 'introduction_cell.png', sep='/'))

# ----------Poster figure 2: Method----------
AP <- list (edge_stroke=1e-5, point_edge_color=alpha ('gray', 0.), fontsize=20, point_fontsize=22/3)
# dataset integration
TB_data$dataset <- gsub ('_[0-9]+$', '', TB_data$paper)
plot_dim_red (TB_data, group.by= c('dataset'), DR='pca' , dims=c(1,2),
                    return_sep=T, nudge_ratio=0.05, plot_type='dim_red_sim',
                    seg_color='black', length_ratio=0.1, AP=AP, nudge_ortho=0.5) [[1]] + 
theme (legend.position='none') +
ggtitle ('Dataset integration') -> plotlist [['method_data']] 

# trajectory construction
meta <- TB_data@meta.data
meta <- meta [!is.na(meta$MGP_PT) & !meta$broad_type %in% c('EPI', 'PE'),]
ggplot () +
plot_two_branch (meta, 
                 select_cells1 = meta$broad_type != 'STB',
                 select_cells2 = meta$broad_type != 'EVT',
                 text_offset=0.005, off_set= 8, thickness=0.01,
                 axis_label='pseudotime', shorten_ratio=0.6, 
                 AP=AP, show_tick=F
) + theme (legend.position='none') +
ggtitle ('Trajectory contruction') -> plotlist[['method_traj']]
plotlist$method_traj

# GRN inference
color_row <- read.csv ( paste ( save_dir4,  'WGCNA/module_genes.csv' , sep='/'), row.names=1)
gene_list <- lapply (as.list (colnames (color_row) ), function (x) {
                             unique (color_row [, x]) })
names (gene_list) <- colnames (color_row)
plot_gene <- list (TB=gene_list$GC1, CTB=gene_list$GC6, STB=gene_list$GC3, EVT=gene_list$GC10)
markers <- find_DE_genes (all_data, save_dir1, group.by='broad_type', label='all_vivo')
data (TF)
tf_data <- all_data [TF, all_data$date!='in_vitro']
pnet <- custom_net_diff_nets (tf_data, plot_gene, markers, nudge_ratio=0.3, size_thres=0.2, AP=AP)
pnet <- lapply (pnet, function (gx){gx+coord_cartesian (clip='off')})
plotlist [['method_GRN']] <- pnet[[1]]+ggtitle ('Trophoblast GRN') + theme (legend.position='none')

# hTSC
highlight <- all_data$assigned_cluster %in% CT$in_vitro_cells
plot_dim_red (all_data, group.by= c('broad_type'), DR='pca' , return_sep=T,
                    size_highlight=highlight, highlight_ratio=4, dims=c(1,2),
                    nudge_ratio=0.05, plot_type='dim_red_sim', seg_color='black', 
                    label_col='none', AP=AP, length_ratio=0.1, nudge_ortho=0.5)[[1]]+
                 theme (legend.position='none') +
                 ggtitle ('In vitro comparison') -> plotlist [['method_hTSC']]

# ----------Poster figure 3: integration----------
plot_dim_red (TB_data, group.by= c('broad_type', 'date'), DR='pca' , dims=c(1,2),
              return_sep=T, nudge_ratio=0.05, plot_type='dim_red_sim',
              seg_color='black', nudge_ortho=0.5, 
              length_ratio=0.1, AP=AP) -> p1

plotlist [['integrate_type']] <- p1[[1]] + labs (fill='')
plotlist [['integrate_date']] <- p1[[2]] + labs (fill='')

data (lineage_markers)
show_genes <- lineage_markers [!names (lineage_markers) %in% c('STR', 'PE')]
show_genes <- show_genes [!show_genes %in% c('PRDM14', 'ELF5', 'ETS2', 'HAND1', 'MMP9', 'KLF17', 'CGB2', 'LHB')]
seurat_heat (TB_data, color_row=show_genes, group.by = c('broad_type'), 
             slot_data='data', heat_name='norm count', center_scale=T,
             column_legend_labels=c('cell type'), row_scale=T, 
             main_width=24, main_height=16, automatic=F,
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
labs (fill='cell type') +theme (legend.position='none', aspect.ratio=0.8) -> plotlist [['traj_gp']]

# pseudotime correlation
metadata <- metadata [!is.na(metadata$MGP_PT),]
plotlist [['traj_pt_cor']] <-pseudo_real_time (metadata, 'date', 'MGP_PT', 'broad_type', lower_b=0, AP=AP) +
        theme (legend.position='bottom')

# signaling pathway
module_scores <- get_module_score (all_data,
                              save_path=paste (sup_save_dir2, 'Data_module_scores.csv', sep='/'))
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

# TF network
plotlist [['traj_network']] <- ggpubr::ggarrange (plotlist=pnet[2:4], nrow=1, 
                                                  common.legend=T, legend='right')

# ----------figure 5: compare ----------
# gene expression pattern
in_vitro <- all_data [, all_data$date == 'in_vitro' ]
plot_genes <- c('POU5F1', 'SOX2', 'NANOG', 'KRT7', 'KRT18', 'TFAP2C',
                'TFAP2A','GATA3','GATA2', 'ELF5', 'CDX2', 'TEAD4', 'HLA-A',
                'HLA-B')
seurat_violin (in_vitro, features=plot_genes, group.by='final_cluster',
               num_col=7, box_plot=F, free_xy='free_y', lower_b=0, AP=AP)+
             theme (aspect.ratio=1.2, legend.position='none')-> plotlist [['hTSC_gene']]

# in GPLVM
plot_met %>% drop_na (hTSC_OKAE) %>% normalize_prob(c('hTSC_OKAE', 'hTSC_TURCO')) %>%
        dim_red_3D_traj ('PT1', 'PT2', 'PT3', c('hTSC_OKAE', 'hTSC_TURCO'), epg, 'x', 'y',
        'z', 'branch_name', all_theta=50, all_phi=0, show_label=T, further_repel=T,
        repel_force=0.2, lab_just=c(0.08, 0.02, 0.02), label_col='broad_type',
        label_traj_text=label_epg, num_col=2, AP=AP, hor_just=0.1, dim_vjust=3)+ 
        labs (fill='relative probability') +
        scale_color_manual (values=rep('black',3))+guides(color=F) -> plotlist [['hTSC_gp']]

# TF activity
module_score <- get_module_score (all_data, append_meta=T, paste (save_dir4, 
                'WGCNA/Data_module_score.csv', sep='/'), pgenes=gene_list)
vitro_mod <- module_score [, module_score$date=='in_vitro']
sel_genes <- rownames (module_score)
names (sel_genes) <- c('TB', 'non-specific', 'STB', 'non-specific', 'EPI', 'CTB', 'ICM', 'ICM',
                       'non-specific', 'EVT', 'cleavage')
row_levels <- partial_relevel (names (sel_genes)) %>% levels()
seurat_heat (vitro_mod, group.by=c('broad_type'),
             row_scale=T, color_row= sel_genes,
             row_legend_labels='WGCNA clusters',
             column_legend_labels='cell type',
             cluster_rows=T, heat_name='module score', center_scale=T,
             automatic=F, left_HA=F, slot_data='counts',
             column_title_fontface='plain',
             column_title_side='bottom', AP=AP, top_HA=F,
             main_width=4, main_height=14, column_rotation=90
) -> plotlist [['hTSC_TF']]

# signaling pathway
sup_save_dir2 <- paste (root, 'manuscript/figureS2', sep='/')
module_scores <- get_module_score (all_data, save_path=paste (sup_save_dir2, 
                                        'Data_module_scores.csv', sep='/'))
colnames (module_scores) <- gsub ('^X', '', colnames (module_scores))
rownames (module_scores) <- gsub ('\\.', '-', rownames (module_scores) )
meta_data <- all_data@meta.data [match (colnames (module_scores), colnames (all_data) ), ]
module_signal <- Seurat::CreateSeuratObject ( module_scores, meta.data = meta_data )
module_signal <- module_signal [, !module_signal$broad_type %in% c('EPI', 'PE', 'hESC', 'hESC-YAN')]
module_signal <- module_signal [, !module_signal$final_cluster %in% c('uCTB')]

sel_seurat <- average_by_group (module_signal, 'final_cluster', rownames (module_signal))
seurat_heat (sel_seurat, 'final_cluster', rownames (module_signal),
             main_width=12, main_height=14, column_split=NA,
             column_rotation=90, show_column_names=T, cluster_column=T,
             center_scale=T, column_legend_labels=c('cell type'), row_scale=T,
             grid_height=5, heat_grid_height=9, automatic=F, AP=AP, top_HA=F,
             heat_name='module score'
)-> plotlist [['hTSC_signal']]

# ----------combine all plots----------
poster_text <- rjson::fromJSON(file = paste (save_dir, 'poster_text.json', sep='/') )
plotlist_t <- c(plotlist, poster_text)
config_arr <- read.csv (paste (save_dir, 'poster_config.csv', sep='/'))
arrange_poster (plotlist_t, paste (save_dir, 'figure4.pdf', sep='/'), config_arr, plot_width=7.5/3, plot_height=7/12)

