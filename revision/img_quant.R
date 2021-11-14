library (tidyverse)
devtools::load_all ('..', export_all=T)
root<- '/mnt/c/Users/Yutong/Documents/bioinformatics/reproduction/results/'
save_dir <- paste (root, 'manuscript/figure4/IF', sep='')
# where the summary data file is located
all_data <- read.csv (paste (root, 'imaging/data/proc_data/all_cleaned2.csv', 
                             sep=''), row.names=1)
all_data %>% filter (experimenter=='Dylan') -> all_data

# ----------Individual pathway comparison----------
all_paths <- c('EGF', 'FGF', 'PD03', 'CHIR', 'FK', 'FK + PD03', 'PD03 + CHIR')

cairo_pdf (paste (save_dir, 'gene_per_condition.pdf', sep='/'), onefile=T)
for (one_path in all_paths){
        print (plot_xy_per_path (all_data, one_path, logxy=F))
        print (boxplot_gene (all_data, one_path, logy=T))
}
dev.off()

# p values
pval_all <- compare_average (all_data , 'condition', ref='base')
write.csv (pval_all, paste (save_dir, 'pval_to_base.csv', sep='/'))

all_data %>% filter (Date!='pd03_siRNA') %>%
        filter (!condition %in% c('PD03 + SB203', 'PD03 + EGF + SB203', 
                                  'PD03 + EGF', 'EGF + SB203', 'PD03 + CHIR')) -> all_data

# fluorescence value
for (one_gene in c('HLAG', 'CGB', 'TFAP2C')){
        cairo_pdf (paste (save_dir, '/Scatter_', one_gene, '_lightgray.pdf', sep=''), width=5, height=7)
        print (fluoro_plot (all_data, one_gene))
        dev.off()
}

# compare absolute cell number
metric<- 'perc'
cairo_pdf (paste (save_dir, '/Box_HLAG_', metric, '.pdf', sep=''), width=5, height=7)
frequency_plot (all_data, 'HLAG', plot_type='error', add_pval='PD03', 
                perc=metric=='perc', sum_table=F)
dev.off()

cairo_pdf (paste (save_dir, '/Box_CGB_', metric, '.pdf', sep=''), width=5, height=7)
frequency_plot (all_data, 'CGB', plot_type='error', add_pval=c('EGF', 'FK', 'PD03'), 
                perc=metric=='perc', sum_table=T)
dev.off()

cairo_pdf (paste (save_dir, '/Box_TFAP2C_', metric, '.pdf', sep=''), width=5, height=7)
frequency_plot (all_data, 'TFAP2C', plot_type='error',  
                perc=metric=='perc', sum_table=T)
dev.off()

# cell number table
table (all_data$condition, all_data$HLAG_pos)

# ----------All pathway comparisons----------
data (format_conf)
pval_all %>% select (feature, mean_expr, condition) %>%
        spread (feature, mean_expr) %>%
        ggplot (aes (x=HLAG, y=CGB)) +
        geom_point ()+
        ggrepel::geom_text_repel (aes (label=condition), size=format_conf$point_fontsize) +
        TBdev::theme_TB ('dotplot', rotation=0)

all_data %>% dplyr::select (dplyr::all_of (c('HLAG', 'CGB', 'TFAP2C', 'condition'))) %>%
        gather ('gene', 'expr', -condition) %>%
        ggplot (aes (x=condition, y=expr)) +
        geom_boxplot (aes (fill=condition))+
        facet_wrap (~gene, scales='free_y') +
        TBdev::theme_TB ('dotplot', rotation=90, AP=AP)+
        theme (panel.border=element_rect (color='grey', fill=NA))+
        scale_y_log10()

# ----------siRNA in differentiating condition----------
pd03 <- all_data %>% filter (Date=='pd03_siRNA') %>% filter (TFAP2C_pos)
metric <- 'perc'
frequency_plot (pd03, 'HLAG', box_plot=T, perc=metric=='perc', sum_table=F,
                add_pval='GATA2', control_group='GFP')
frequency_plot (pd03, 'CGB', box_plot=T, perc=metric=='perc', sum_table=F,
                add_pval='GATA2', control_group='GFP')

frequency_plot (all_data, 'CGB', box_plot=T, perc=metric=='perc', sum_table=F,
                add_pval='PD03', control_group='GFP')+theme (aspect.ratio=0.7)

# ----------PCA----------
select_cols <- c('HLAG', 'TFAP2C', 'CGB', 'area', 'major_axis_length',
                 'minor_axis_length', 'nuc_circ', 'perimeter')
seurat_ob <- Seurat::CreateSeuratObject (t(all_data [, select_cols]), 
                                         meta.data=all_data)
seurat_ob[['RNA']]@data <-seurat_ob[['RNA']]@counts
Seurat::VariableFeatures (seurat_ob) <- rownames (seurat_ob)
seurat_ob <- Seurat::ScaleData(seurat_ob)
seurat_ob <- Seurat::RunPCA (seurat_ob, npc=nrow(seurat_ob))

seurat_ob$logHLAG <- log (seurat_ob$HLAG+1)
seurat_ob$logCGB <- log (seurat_ob$CGB+1)
seurat_ob$genes <- ifelse (seurat_ob$HLAG_pos, 'HLA-G', 'none')
seurat_ob$genes [seurat_ob$CGB_pos] <- 'CGB'
seurat_ob$genes [seurat_ob$CGB_pos & seurat_ob$HLAG_pos] <- 'CGB & HLA-G'

TBdev::plot_dim_red (seurat_ob, group.by=c('condition', 'genes'), further_repel=F, dims=1:3, all_theta=-50, all_phi=-30)
TBdev::plot_dim_red (seurat_ob, group.by='genes', further_repel=F, dims=1:3, all_theta=-50, all_phi=-30)
TBdev::plot_gene_PC (seurat_ob)
