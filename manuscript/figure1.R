# generate the figure 1 of the manuscript
# unified atlas of trophoblast development

setwd ('..')
devtools::load_all ()
library (tidyverse)
library (GOSemSim)

root_dir <- '/mnt/c/Users/Yutong/Documents/bioinformatics/reproduction/'
root <- paste (root_dir, 'results/', sep='/')
merge_dir <- paste (root, 'XLYBPZ_Dylan_dir', sep='/')
save_dir <- paste (root, 'manuscript/figure1', sep='/')
sup_save_dir <- paste (root, 'manuscript/figureS1', sep='/')
x <- load (paste (merge_dir, 'final_merged_vivo.Robj', sep='/') )
all_data <- get (x)

data (CT)
in_vivo <- all_data [, all_data$date != 'in_vitro' & !(all_data$broad_type %in% CT$in_vitro_cells)]
tb_only <- in_vivo [, !in_vivo $broad_type %in% CT$non_TB_lineage]

# key genes
plot_genes <- c('POU5F1', 'SOX2', 'NANOG', 'GATA3', 'CDX2', 'TFAP2C')
# heatmap
markers <- find_DE_genes (in_vivo, save_dir, feature='broad_type', label='all_vivo')
DE_genes <- unique_DE_genes (markers, 6)
DE_genes %>% filter (!group %in% CT$non_TB_lineage) %>% 
        select (group, feature) %>% deframe () -> show_genes

# run GO and KEGG
d <- godata(org.Hs.eg.db, ont="BP")
markers %>% filter (!group %in% CT$non_TB_lineage) -> tb_markers
kk <- compare_cluster_enrichment (tb_markers, d, enrich_area='KEGG')

p1 <- plot_dim_red (in_vivo, group.by= c('broad_type'), DR='pca', all_labels=T,
                    return_sep=F, nudge_ratio=0.2)
p2 <- seurat_violin (in_vivo, features=plot_genes, group.by='broad_type')
p3 <- seurat_heat (tb_only, color_row=show_genes, group.by = c('broad_type', 'date'), 
                       slot_data='data', heat_name='norm count', center_scale=T,
                       column_legend_labels=c('cell type', 'date'), row_scale=T)

set.seed (100)
p4 <- display_cluster_enrichment (kk, show_graph='emap', feature_vec= tb_markers$group, 
                                     show_num=20, vert_just=2.2) + labs (fill='')
psea_df <- run_GSEA_all_types (tb_markers, org.Hs.eg.db, save_path = 
                                  paste (sup_save_dir, 'PSEA_broad_type_all_vivo.csv', sep='/'))
set.seed (100)
p5 <- ridge_all_types (psea_df, show_num=40, not_show_prop=0.01) + labs (fill='') 

grob_list <- list (p1, p2 + xlab ('') + ylab ('norm count')+
                   theme (axis.title.x=element_blank (), aspect.ratio=1.5) , 
                   p3, p4, p5+theme (aspect.ratio=2.5))
lay_mat <- cbind (c(1,3,4), c(2,3,5))
arrange_plots (grob_list, paste (save_dir, 'final_figure1.pdf', sep='/'), lay_mat )
