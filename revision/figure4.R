# TF network in vivo and in vitro
# clonogenicity assay to verify the network

devtools::load_all ('..', export_all = F)
library (tidyverse)

root_dir <- '/mnt/c/Users/Yutong/Documents/bioinformatics/reproduction/'
root <- paste (root_dir, 'results/', sep='/')
merge_dir <- paste (root, 'XLYBPZ_Dylan_dir', sep='/')
#save_dir <- paste (root, 'manuscript/figure2', sep='/')
#sup_save_dir <- paste (root, 'manuscript/figureS2', sep='/')
save_dir <- paste (root, 'manuscript/figure4', sep='/')
all_data <- readRDS(paste (merge_dir, 'merged_blastoid.rds', sep='/') )

# ----------aesthetics----------
data (format_conf)
new_order <- c(format_conf$cell_order, 'hESC', 'nESC', 'bEPI-YANA', 'hTSC-OKAE','bTSC', 'bTB-YANA')
blastoid_color <- c('nESC'='#71d300', 'bTB-YANA'='#ff1900', 'bTSC'='#ff8500', 'bEPI-YANA'='#00ff9d')
new_color <- list (color_vec=c(format_conf$color_vec, blastoid_color), cell_order=new_order)
new_color$color_vec <- new_color$color_vec [order (partial_relevel (names (
                                        new_color$color_vec), new_color$cell_order))]

# ----------Figure 4A: WGCNA nets----------
color_row <- read.csv ( paste ( save_dir,  'WGCNA/module_genes.csv' , sep='/'), row.names=1)
gene_list <- lapply (as.list (colnames (color_row) ), function (x) {
                             unique (color_row [, x]) })
names (gene_list) <- colnames (color_row)

# load data
data (TF) #load a vector of TF gene names
TF_WG <- all_data[TF,all_data$date != 'in_vitro']
fil_vivo <- filter_genes (TF_WG, 0.2)
rm (TF_WG)

save_dir1 <- paste (root, 'manuscript/figure1', sep='/')
markers <- find_DE_genes (all_data, save_dir1, group.by='broad_type', label='all_vivo')
plot_gene <- list (TB=gene_list$GC1, CTB=gene_list$GC6, STB=gene_list$GC3, EVT=gene_list$GC10)

plotlist <- custom_net_diff_nets (fil_vivo, plot_gene, markers, nudge_ratio=0.3, size_thres=0.2)
plotlist2 <- lapply (plotlist, function(gx){gx+theme(aspect.ratio=0.85)+coord_cartesian (clip='off')})
p1 <- ggpubr::ggarrange (plotlist=plotlist2, nrow=1, common.legend=T, legend='right')
p1

# ----------figure 4B: in vivo heatmap----------
module_score <- get_module_score (all_data, append_meta=T, paste (save_dir, 
                'WGCNA/Data_module_score.csv', sep='/'), pgenes=gene_list)
vivo_mod <- module_score [, module_score$date != 'in_vitro']
vivo_scale <- scale_seurat (vivo_mod, row_scale=T, slot_data='counts')
sel_genes <- rownames (vivo_scale)
names (sel_genes) <- c('TB', 'non-specific', 'STB', 'non-specific', 'EPI', 'CTB', 'ICM', 'ICM',
                       'non-specific', 'EVT', 'cleavage')
p2 <- seurat_heat (vivo_scale, group.by=c('broad_type', 'date'),
                 row_scale=T, color_row= sel_genes,
                 column_reorder_levels = list (format_conf$branch_order, format_conf$cell_order),
                 column_legend_labels= c('cell type', 'date'), 
                 row_legend_labels='WGCNA clusters',
                 cluster_rows=T, heat_name='module score', center_scale=T,
                 column_rotation=90,
                 main_width=21, main_height=14,
                 automatic=F, left_HA=F, slot_data='counts'
)

# ----------figure 4C: in vitro heatmap----------
vitro_mod <- module_score [, module_score$broad_type %in% c('hESC', 'hESC-YAN', 
                                        'bEPI-YANA', 'hTSC-OKAE', 'bTB-YANA')]
sel_genes_vitro <- c('GC11', 'GC5', 'GC7', 'GC8', 'GC1', 'GC10', 'GC3', 'GC6',
                     'GC2', 'GC4', 'GC9')
names (sel_genes_vitro) <- names (sel_genes) [match (sel_genes_vitro, sel_genes)]
p3 <- seurat_heat (vitro_mod, group.by=c('broad_type'),
                 row_scale=T, color_row= sel_genes_vitro,
                 row_legend_labels='WGCNA clusters',
                 column_legend_labels='cell type',
                 cluster_rows=F, heat_name='norm count', center_scale=T,
                 automatic=F, left_HA=F, slot_data='counts',
                 column_title_fontface='plain',
                 column_title_side='top',
                 main_width=6, main_height=14, column_rotation=90, AP=new_color
)
p3

# ----------Hub genes----------
save_dir3 <- paste (root, 'manuscript/figure3', sep='/')
markers <- find_DE_genes (all_data [, all_data$date=='in_vitro'], 
                          save_dir3, group.by='broad_type', 
                          label='pairwise_blastoid', method='pairwise')
# EVT genes
tf_markers <- markers [markers$feature %in% gene_list$GC10,]
evt_hub <- read.csv (paste (save_dir, 'WGCNA/TF_EVT.csv',sep='/'))
tf_markers$hub_score <- evt_hub$kWithin [match (tf_markers$feature, evt_hub$feature)]
label_genes <- evt_hub$feature [base::order (evt_hub$kWithin, decreasing=T)][1:20]
seurat_volcano (tf_markers, 'bTB-YANA', 'hTSC-OKAE', length_ratio=0.5, label_genes=label_genes,
                pval='hub_score', logy=F, show_gene_num=14)+ylab ('EVT hub score')-> p4

# CTB genes
tf_markers <- markers [markers$feature %in% gene_list$GC6,]
ctb_hub <- read.csv (paste (save_dir, 'WGCNA/TF_CTB.csv',sep='/'))
ctb_hub <- ctb_hub [ctb_hub$cluster=='GC6',]
tf_markers$hub_score <- ctb_hub$kWithin [match (tf_markers$feature, ctb_hub$feature)]
label_genes <- ctb_hub$feature [base::order (ctb_hub$kWithin, decreasing=T)][1:20]
seurat_volcano (tf_markers, 'bTB-YANA', 'hTSC-OKAE', length_ratio=0.5, label_genes=label_genes,
                pval='hub_score', logy=F)+ylab ('CTB hub score')-> p5

# ----------arrange all----------
grob_list <- list (p1, p2, p3, p4+theme (aspect.ratio=0.7), 
                   p5+theme (aspect.ratio=0.7), ggplot ()+theme_minimal())
lay_mat <- matrix(c(1,1,1,
                    2,2,3,
                    2,2,3,
                    4,6,6,
                    5,6,6
                    ),
                  nrow=3) %>% t()
arrange_plots (grob_list, paste (save_dir, 'figure4_upper.pdf', sep='/'),
               lay_mat, plot_width=7, plot_height=4)

# ----------net for siRNA----------
library (igraph)
# siRNA KD counts
read.csv (paste (save_dir, 'clonogenicity/clonogenicity_KD_efficacy.csv', sep='/')) %>% 
        select (!X) %>% deframe () -> gene_expr 
gene_expr <- 1 - gene_expr # measure by the degree of downregulation

# plot network
plot_gene <- c('CTNNB1', 'GATA2', 'GATA3', 'GCM1', 'GRHL1', 'HMGA1', 'HMOX1',
               'IRX4', 'MAZ', 'MSX2', 'NFE2L3', 'NR2F2', 'OVOL1', 'PCBP2',
               'PITX1', 'SP6', 'SSRP1', 'TEAD3', 'TFAP2C', 'TFEB', 'TRIM28',
               'ZFHX3', 'ZNF750')
graph_net <- plot_WGCNA_net (fil_vivo, plot_gene, return_igraph=T, thres=0., scale_node_size=0.0001)
#graph_net <- add_pseudotime_to_net (fil_vivo, graph_net)
V(graph_net)$size <- gene_expr [match (attr (V(graph_net), 'names'), names (gene_expr) )]*200
limits <- range (fil_vivo [plot_gene,][['RNA']]@data)

cairo_pdf (paste (save_dir, 'siRNA_net.pdf', sep='/'), height=4, width=6)
custom_net (graph_net, coloring=NULL, 
            size_legend=guide_legend(),
            order_leftright=NULL, plot_title='siRNA') +coord_cartesian (clip='off')
dev.off ()

in_vitro <- all_data [, all_data$date=='in_vitro']
cairo_pdf (paste (save_dir, 'vivo_siRNA.pdf', sep='/'), height=8, width=14)
seurat_heat (fil_vivo, group.by=c('broad_type', 'date'),
                 row_scale=T, color_row= plot_gene, cluster_rows=T,
                 column_legend_labels= c('cell type', 'date'), 
                 heat_name='norm count', center_scale=T,
                 column_rotation=90,
                 main_width=21, main_height=14,
                 automatic=F, left_HA=F, slot_data='counts'
)
dev.off ()

cluster_level <- c('GFP', 'PITX1', 'CTNNB1', 'TFAP2C', 'HMOX1', 'ZNF750',
                   'SP6', 'MSX2', 'OVOL1', 'GRHL1', 'GCM1', 'TEAD3', 'GATA2',
                   'ZFHX3', 'GATA3', 'IRX4', 'NFE2L3', 'NR2F2', 'TFEB',
                   'TRIM28', 'MAZ', 'SSRP1', 'PCBP2', 'HMGA1')
cairo_pdf (paste (save_dir, 'vitro_siRNA.pdf', sep='/'), height=8, width=9)
seurat_heat (in_vitro, group.by=c('broad_type'),
                 row_scale=T, color_row= cluster_level,
                 column_legend_labels= c('cell type'), 
                 heat_name='norm count', center_scale=T,
                 column_rotation=90,
                 main_width=14, main_height=14,
                 automatic=F, left_HA=F, slot_data='counts'
)
dev.off ()

