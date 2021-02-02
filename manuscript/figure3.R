devtools::load_all('..', export_all=F)
library (tidyverse)

root_dir <- '/mnt/c/Users/Yutong/Documents/bioinformatics/reproduction/'
root <- paste (root_dir, 'results/', sep='/')
merge_dir <- paste (root, 'XLYBPZ_Dylan_dir', sep='/')
save_dir <- paste (root, 'manuscript/figure3', sep='/')
sup_save_dir <- paste (root, 'manuscript/figureS3', sep='/')
all_data <- get (load (paste (merge_dir, 'final_merged_tb.Robj', sep='/') ))
all_data <- all_data [, all_data$assigned_cluster != 'uCTB']

# ----------figure 3A ----------
TSC_cells <- all_data$revised %in% c('hTSC_OKAE', 'hTSC_TURCO') 
#all_data$new_cluster2 [TSC_cells] <- as.character (all_data$revised [TSC_cells])
#all_data$new_cluster2 <- ML$partial_relevel (all_data$new_cluster2, ML$cell_order)

data (CT)
highlight <- all_data$assigned_cluster %in% CT$in_vitro_cells
p1 <- plot_dim_red (all_data, group.by= c('assigned_cluster'), DR='pca' , return_sep=F,
                    size_highlight=highlight, highlight_font=2, dims=c(1,2), nudge_ratio=0.2)
p1
# ----------figure 3B ----------
in_vitro <- all_data [, all_data$date == 'in_vitro' ]
data (lineage_markers)
plot_genes <- lineage_markers [names (lineage_markers)  == 'TB']
plot_genes <- c('POU5F1', 'SOX2', 'NANOG', plot_genes, 'TP63')
p2 <- seurat_violin (in_vitro, features=plot_genes, group.by='revised', num_col=5, box_plot=F)

# ----------figure 3C: cell-cell correlation ----------
select_cells2 <- all_data$revised %in% CT$in_vitro_cells
select_cells1 <- !(all_data$revised %in% CT$in_vitro_cells)
all_cor <- compute_all_cor (all_data, method= 'correlation', assay='RNA',
                            select_cells2=select_cells2,
                            select_cells1=select_cells1)
devtools::load_all('..', export_all=F)
p4 <- cell_violin (all_cor, all_data@meta.data, c('assigned_cluster', 'revised'), 
                   box_plot=T, num_col=2, legend_col=1)

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

show_meta %>% slice_min (PT3, n=nrow(show_meta)-3) -> plot_met
p3 <- dim_red_3D_traj (plot_met, 'PT1', 'PT2', 'PT3', c('hTSC_OKAE', 'hTSC_TURCO'), epg, 'x', 'y',
        'z', 'branch_name', all_theta=50, all_phi=0, show_label=T, further_repel=T,
        repel_force=0.2, lab_just=c(0.08, 0.02, 0.02), label_col='broad_type',
        label_traj_text=label_epg) + labs (fill='probability') +
        scale_color_manual (values=rep('black',3))+guides(color=F)

# ----------figure 3E: DE gene expression----------
TSC_data <- all_data [, TSC_cells]
TSC_data$assigned_cluster <- TSC_data$revised
TSC_data$select <- 'select'
all_data$select <- 'non_select'
heat_data <- merge (all_data, TSC_data)
select_cells <- heat_data$select == 'select'

save_dir1 <- paste (root, 'manuscript/figure1', sep='/')
markers <- find_DE_genes ( all_data,  save_dir1, group.by='broad_type', label='all_vivo')
DE_genes <- unique_DE_genes (markers, 5) %>% as.data.frame ()
DE_genes$group <- as.character (DE_genes$group)
DE_genes %>% select (group, feature) %>% deframe () -> show_genes

seurat_param <- list (
        heat_name=c('norm count'),
        column_legend_labels=c('cell type'),
        main_width = c(10, 7),
        column_rotation=c(90, 0),
        show_column_anna=c(T, F),
        column_title_fontface=c('plain','bold'),
        grid_height=5
)
p5 <- seurat_heat_highlight (heat_data, select_cells, show_genes,
                             c('assigned_cluster'), average=T, return_sep=T,
                             seurat_heat_params=seurat_param, row_scale=T
)

p5_final <- ComplexHeatmap::draw (p5[[1]]+p5[[2]], 
                                  ht_gap = unit(c(0.3, 0.5), "cm"))

# ----------figure 3F: TF network----------
invivo <- all_data [, all_data$date != 'in_vitro']
data (TF) #load a vector of TF gene names
TF_WG <- invivo [TF,]
fil_vivo <- filter_genes (TF_WG, 0.2)

save_dir4 <- paste (root, 'manuscript/figure4', sep='/')
color_row <- read.csv ( paste ( save_dir4,  'WGCNA/module_genes.csv' , sep='/'), row.names=1)
gene_list <- lapply (as.list (colnames (color_row) ), function (x) {
                             unique (color_row [, x]) })
names (gene_list) <- colnames (color_row)

merged4_6 <- c(as.character (gene_list$GC6), as.character (gene_list$GC4))
merged4_6 <- merged4_6 [!grepl ('^HNRNP', merged4_6)]

graph_net <- plot_WGCNA_net (fil_vivo, merged4_6, return_igraph=T, thres=0.6)
graph_net <- add_pseudotime_to_net (fil_vivo, graph_net)

save_dir1 <- paste (root, 'manuscript/figure1', sep='/')
markers <- find_DE_genes (all_data, save_dir1, group.by='broad_type', label='all')
celltypes <- c('hTSC_OKAE', 'hTSC_TURCO', 'hESC')
graph_list <- custom_net_cells (markers, graph_net, celltypes,
                                normalize_cell='TB', return_sep=T,
                                hide_node_thres=0.2)

# arrange all figures
grob_list <- c(list (p1 +labs (fill=''), 
                   p2 + ylab('norm count')+ theme (aspect.ratio=1., axis.title.x=element_blank()) +
                           labs (fill =''), p3, 
                   p4+ labs (fill='') + theme (axis.title.x=element_blank(), aspect.ratio=0.5), 
                   p5_final), graph_list)
lay_mat <- matrix(c(1, 1, 1, 2, 2, 2, 
                    3, 3, 3, 3, 6, 6,
                    4 ,4 ,4 ,4 ,7 ,7,
                    5 ,5 ,5 ,5 ,8 ,8 
                    ),
                  nrow=6) %>% t()
arrange_plots (grob_list, paste (save_dir, 'final_figure3.pdf', sep='/'),
               lay_mat, plot_width=3.5, plot_height=8)
