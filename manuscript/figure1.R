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
all_data <- get(load (paste (merge_dir, 'final_merged_tb.Robj', sep='/') ))

data (CT)
#all_data <- all_data [, all_data$assigned_cluster != 'uCTB']
in_vivo <- all_data [, all_data$date != 'in_vitro' & !(all_data$broad_type %in% CT$in_vitro_cells)]
tb_only <- in_vivo [, !in_vivo $broad_type %in% CT$non_TB_lineage]

epi <- in_vivo [, in_vivo$broad_type %in% c('ICM', 'EPI')]
AP <- list (color_vec = c('EPI1'='red', 'EPI2'='green'))
plot_dim_red (epi, group.by= c('assigned_cluster'), DR='pca', return_sep=T,
              nudge_ratio=0.2, plot_type='dim_red', AP=AP)

# key genes
plot_genes <- c('GATA3', 'CDX2', 'TFAP2C', 'SDC1', 'CGB1', 'HLA-G', 'POU5F1', 'NANOG', 'SOX17')
# heatmap
markers <- find_DE_genes (in_vivo, save_dir, group.by='broad_type', label='all_vivo')
DE_genes <- unique_DE_genes (markers, 6, iterative=5)
DE_genes %>% filter (!group %in% CT$non_TB_lineage) %>% 
        dplyr::select (group, feature) %>% deframe () -> show_genes

p1 <- plot_dim_red (in_vivo, group.by= c('broad_type'), DR='pca', return_sep=T,
                    nudge_ratio=0.2, plot_type='dim_red_sim')
p2 <- seurat_violin (in_vivo, features=plot_genes, group.by='broad_type')
p3 <- seurat_heat (tb_only, color_row=show_genes, group.by = c('broad_type', 'date'), 
                       slot_data='data', heat_name='norm count', center_scale=T,
                       column_legend_labels=c('cell type', 'date'), row_scale=T)

# run KEGG
markers <- find_DE_genes (in_vivo, save_dir, group.by='broad_type', label='pairwise', method='pairwise')

tb_ctb <- markers %>% filter (group == 'CTB' & compare_group == 'TB') 
psea_tb_ctb <- run_GSEA_all_types (tb_ctb, org.Hs.eg.db, save_path=paste (save_dir, 'PSEA_TB_CTB.csv', sep='/'))
devtools::load_all ('..', export_all=F)
p4 <- enrich_bar (psea_tb_ctb, org.Hs.eg.db, show_num=5, markers=tb_ctb,
            show_gene_labels=6, extend_axis_pos=1.2, extend_axis_neg=1.2, nudge_x=0.2, shrink_ratio=0.7)

ctb_stb <- markers %>% filter (group == 'STB' & compare_group == 'CTB') 
psea_ctb_stb <- run_GSEA_all_types (ctb_stb, org.Hs.eg.db, save_path=paste (save_dir, 'PSEA_CTB_STB.csv', sep='/'))
p5 <- enrich_bar (psea_ctb_stb, org.Hs.eg.db, show_num=5, markers=ctb_stb,
            show_gene_labels=6, extend_axis_pos=1.2, extend_axis_neg=1.2, nudge_x=0.12, shrink_ratio=0.7)

ctb_evt <- markers %>% filter (group == 'EVT' & compare_group == 'CTB') 
psea_ctb_evt <- run_GSEA_all_types (ctb_evt, org.Hs.eg.db, save_path=paste (save_dir, 'PSEA_CTB_EVT.csv', sep='/'))
p6 <- enrich_bar (psea_ctb_evt, org.Hs.eg.db, show_num=5, markers=ctb_evt,
            show_gene_labels=6, extend_axis_pos=1.2, extend_axis_neg=3, nudge_x=0.2, shrink_ratio=0.7)

grob_list <- list (p1[[1]]+labs (fill=''), p2 + xlab ('') + ylab ('norm count')+
                   theme (axis.title.x=element_blank (), aspect.ratio=1) , 
                   p3, p4 +theme (aspect.ratio=1.5), 
                   p5+theme (aspect.ratio=1.5), p6+theme (aspect.ratio=1.5))
lay_mat <- matrix(c(1, 1, 1, 2, 2, 2, 
                    3, 3, 3, 3, 3, 3,
                    4, 4, 5, 5, 6, 6),
                  nrow=6) %>% t()
arrange_plots (grob_list, paste (save_dir, 'final_figure1.pdf', sep='/'), lay_mat, plot_width=3)

# save individual plots
#all_plots <- list (ctb=p4, stb=p5, evt=p6)
#
#for (i in 1:length(all_plots)){
#        cairo_pdf (paste (save_dir, '/GSEA_', names (all_plots)[i], '.pdf',sep=''), width=8, height=8)
#        print (all_plots[[i]])
#        dev.off ()
#}

