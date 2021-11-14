# generate the figure 1 of the manuscript
# unified atlas of trophoblast development

devtools::load_all ('..', export_all=F)
#library (TBdev)
library (tidyverse)
library (org.Hs.eg.db)

root_dir <- '/mnt/c/Users/Yutong/Documents/bioinformatics/reproduction/'
root <- paste (root_dir, 'results/', sep='/')
merge_dir <- paste (root, 'XLYBPZ_Dylan_dir', sep='/')
save_dir <- paste (root, 'manuscript/figure1', sep='/')
sup_save_dir <- paste (root, 'manuscript/figureS1', sep='/')
all_data <- readRDS(paste (merge_dir, 'revision_all.rds', sep='/') )

data (CT)
#all_data <- all_data [, all_data$assigned_cluster != 'uCTB']
in_vivo <- all_data [, all_data$date != 'in_vitro' & !(all_data$broad_type %in% CT$in_vitro_cells)]
tb_only <- in_vivo [, !in_vivo $broad_type %in% CT$non_TB_lineage]

# key genes
plot_genes <- c('GATA3', 'CDX2', 'TFAP2C', 'SDC1', 'CGB1', 'HLA-G', 'POU5F1', 'NANOG', 'SOX17')
# heatmap
markers <- find_DE_genes (in_vivo, save_dir, group.by='broad_type', label='all_vivo')
DE_genes <- unique_DE_genes (markers, 6, iterative=5)
DE_genes %>% filter (!group %in% CT$non_TB_lineage) %>% 
        dplyr::select (group, feature) %>% deframe () -> show_genes

p1 <- plot_dim_red (in_vivo, group.by= c('broad_type'), DR='pca', return_sep=T,
                    nudge_ratio=0.2, plot_type='dim_red_sim', nudge_ortho=0.7)
p2 <- seurat_violin (in_vivo, features=plot_genes, group.by='broad_type', lower_b=0)
p3 <- seurat_heat (tb_only, color_row=show_genes, group.by = c('broad_type', 'date'), 
                       slot_data='data', heat_name='norm count', center_scale=T,
                       column_legend_labels=c('cell type', 'date'), row_scale=T)

# run KEGG
markers <- find_DE_genes (in_vivo, save_dir, group.by='broad_type', label='pairwise', method='pairwise')
psea_tb_ctb <- run_GSEA_pairwise (markers, org.Hs.eg.db, paste (save_dir, 'PSEA_TB_CTB.csv', sep='/'),
                   group1='CTB', group2='TB')
rich_forest (psea_tb_ctb, markers, org.Hs.eg.db, show_num=5, show_gene_labels=6,
             shrink_ratio=0.7, nudge_x=0.2, extend_x_pos=1.2, extend_x_neg=1.2) -> p4

psea_ctb_stb <- run_GSEA_pairwise (markers, org.Hs.eg.db, paste (save_dir, 'PSEA_CTB_STB.csv', sep='/'),
                   group1='STB', group2='CTB')
rich_forest (psea_ctb_stb, markers, org.Hs.eg.db, show_num=5, show_gene_labels=6,
             shrink_ratio=0.7, nudge_x=0.2, extend_x_pos=1.2, extend_x_neg=1.2) -> p5

psea_ctb_evt <- run_GSEA_pairwise (markers, org.Hs.eg.db, paste (save_dir, 'PSEA_CTB_EVT.csv', sep='/'),
                   group1='EVT', group2='CTB')
rich_forest (psea_ctb_evt, markers, org.Hs.eg.db, show_num=5, show_gene_labels=6,
             shrink_ratio=0.7, nudge_x=0.2, extend_x_pos=1.2, extend_x_neg=1.2) -> p6

grob_list <- list (p1[[1]]+labs (fill=''), p2 + xlab ('') + ylab ('norm count')+
                   theme (axis.title.x=element_blank (), aspect.ratio=1) , 
                   p3, p4 +theme (aspect.ratio=1.5), 
                   p5+theme (aspect.ratio=1.5), p6+theme (aspect.ratio=1.5))
lay_mat <- matrix(c(1, 1, 1, 2, 2, 2, 
                    3, 3, 3, 3, 3, 3,
                    4, 4, 5, 5, 6, 6),
                  nrow=6) %>% t()
arrange_plots (grob_list, paste (save_dir, 'final_figure1.pdf', sep='/'), lay_mat, plot_width=3)
save_indiv_plots (grob_list, paste (save_dir, 'figure1', sep='/'), lay_mat, plot_width=3)

# save csv
write.csv (as.matrix (all_data [, all_data$date!='in_vitro'] [['RNA']]@data), 
           paste (merge_dir, 'final_merged.csv', sep='/') )

data (TF, package='TBdev')
write.csv (TF, paste (save_dir, 'TF_list.csv', sep='/') )
