setwd ('..')
devtools::load_all()
library (tidyverse)

root_dir <- '/mnt/c/Users/Yutong/Documents/bioinformatics/reproduction/'
root <- paste (root_dir, 'results/', sep='/')
merge_dir <- paste (root, 'XLYBPZ_Dylan_dir', sep='/')
save_dir <- paste (root, 'manuscript/figure3', sep='/')
sup_save_dir <- paste (root, 'manuscript/figureS3', sep='/')
x <- load (paste (merge_dir, 'final_merged_tb.Robj', sep='/') )
all_data <- get (x)

# ----------figure 3A ----------
TSC_cells <- all_data$revised %in% c('hTSC_OKAE', 'hTSC_TURCO') 
#all_data$new_cluster2 [TSC_cells] <- as.character (all_data$revised [TSC_cells])
#all_data$new_cluster2 <- ML$partial_relevel (all_data$new_cluster2, ML$cell_order)

data (CT)
highlight <- all_data$assigned_cluster %in% CT$in_vitro_cells
p1 <- plot_dim_red (all_data, by_group = c('assigned_cluster'), DR='pca' , return_sep=F,
                    size_highlight=highlight, highlight_font=2, dims=c(1,2), nudge_ratio=0.2)
p1
# ----------figure 3B ----------
in_vitro <- all_data [, all_data$date == 'in_vitro' ]
data (lineage_markers)
plot_genes <- lineage_markers [names (lineage_markers)  == 'TB']
plot_genes <- c('POU5F1', 'SOX2', 'NANOG', plot_genes, 'TP63')
p2 <- seurat_violin (in_vitro, features=plot_genes, group.by='revised', num_col=5)

# ----------figure 3C----------
sup_save_dir2 <- paste (root, 'manuscript/figureS2', sep='/')
epg <- read.csv (paste (sup_save_dir2, 'result/STREAM_graph.csv', sep='/'))
show_meta <- all_data@meta.data
show_meta %>% filter (!is.na (MGP_PT) ) %>% filter (!broad_type %in% c('PE', 'EPI') ) -> show_meta
p3 <- dim_red_3D_traj (show_meta, 'PT1', 'PT2', 'PT3', c('hTSC_OKAE', 'hTSC_TURCO'), epg, 'x', 'y',
        'z', 'branch', all_theta=50, all_phi=0, show_label=T, further_repel=T,
        repel_force=0.2, lab_just=c(0.08, 0.02, 0.02), label_col='broad_type') + 
        labs (fill='probability')

# ----------figure 3D: cell-cell correlation ----------
select_cells2 <- all_data$revised %in% CT$in_vitro_cells
select_cells1 <- !(all_data$revised %in% CT$in_vitro_cells)
all_cor <- compute_all_cor (all_data, method= 'correlation', assay='RNA',
                            select_cells2=select_cells2,
                            select_cells1=select_cells1)
p4 <- cell_violin (all_cor, all_data@meta.data, c('assigned_cluster', 'revised'))

# ----------figure 3E: DE gene expression----------
TSC_data <- all_data [, TSC_cells]
TSC_data$assigned_cluster <- TSC_data$revised
TSC_data$select <- 'select'
all_data$select <- 'non_select'
heat_data <- merge (all_data, TSC_data)
select_cells <- heat_data$select == 'select'

save_dir1 <- paste (root, 'manuscript/figure1', sep='/')
markers <- find_DE_genes ( all_data,  save_dir1, feature='broad_type', label='all_vivo')
DE_genes <- unique_DE_genes (markers, 5) %>% as.data.frame ()
DE_genes$group <- as.character (DE_genes$group)
DE_genes %>% select (group, feature) %>% deframe () -> show_genes

seurat_param <- list (
        heat_name=c('norm count'),
        column_legend_labels=c('cell type'),
        main_width = c(10),
        column_rotation=c(90, 0),
        show_column_anna=c(T, F),
        grid_height=5
)
p5 <- seurat_heat_highlight (heat_data, select_cells, show_genes,
                             c('assigned_cluster'), average=T, return_sep=T,
                             seurat_heat_params=seurat_param, row_scale=T)

p5_final <- ComplexHeatmap::draw (p5[[1]]+p5[[2]], 
                                  ht_gap = unit(c(0.3, 0.5), "cm"))

# arrange all figures
grob_list <- list (p1 +labs (fill=''), 
                   p2 + ylab('norm count')+ theme (aspect.ratio=1., axis.title.x=element_blank()) +
                           labs (fill =''), p3, 
                   p4+ labs (fill='') + theme (axis.title.x=element_blank()), 
                   p5_final)
lay_mat <- matrix(c(1, 1, 1, 2, 2, 2, 
                    3, 3, 3, 3, 4, 4,
                    5 ,5 ,5 ,5 ,4 ,4 
                    ),
                  nrow=6) %>% t()
arrange_plots (grob_list, paste (save_dir, 'final_figure3.pdf', sep='/'),
               lay_mat, plot_width=3.5, plot_height=8)
