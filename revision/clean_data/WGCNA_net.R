# Identify the significant transcription factor (TF) network for siRNA screenining
root_dir <- '/mnt/c/Users/Yutong/Documents/bioinformatics/reproduction/'
root <- paste (root_dir, 'results/', sep='/')
merge_dir <- paste (root, 'XLYBPZ_Dylan_dir', sep='/')
save_dir <- paste (root, 'manuscript/figure4', sep='/')
sup_save_dir <- paste (root, 'manuscript/figureS4', sep='/')

setwd ('..')
devtools::load_all ('.', export_all=T)
library (WGCNA)

all_data <- get(load (paste (merge_dir, 'final_merged_tb.Robj', sep='/') ))
# construct TF network using the in vivo set only
invivo <- all_data [, all_data$date != 'in_vitro']
data (TF) #load a vector of TF gene names
TF_WG <- invivo [TF,]
fil_vivo <- filter_genes (TF_WG, 0.2)

# ----------run WGCNA----------
find_eigengene (fil_vivo, save_dir, minModuleSize =35, cluster_num='all',
                top_markers='all', return_eigengene=F)

color_row <- read.csv ( paste ( save_dir,  'WGCNA/module_genes.csv' , sep='/'), row.names=1)
gene_list <- lapply (as.list (colnames (color_row) ), function (x) {
                             unique (color_row [, x]) })
names (gene_list) <- colnames (color_row)
lapply (gene_list, length)
plot_all_WGCNA_nets (fil_vivo, gene_list, paste (save_dir, 'WGCNA/GC_maps.pdf', sep='/'))

# ----------module scores----------
module_score <- get_module_score (fil_vivo, paste (save_dir, 
                'WGCNA/Data_module_score.csv', sep='/'), pgenes=gene_list)
colnames (module_score) <- gsub ('^X', '', colnames (module_score))

show_meta <- all_data@meta.data [!is.na (all_data$MGP_PT) & !all_data$broad_type %in% c('EPI', 'PE') ,]
module_score <- module_score [, match (rownames (show_meta), colnames (module_score) ) ]
module_data <- Seurat::CreateSeuratObject ( module_score, meta.data = show_meta)
module_data <- module_data [, !is.na (module_data$MGP_PT) ]
sel_genes <- rownames (module_data)

data (format_conf)
devtools::load_all ('.', export_all=F)
seurat_heat (module_data, group.by=c('epil_branch','broad_type'),
                 row_scale=T, color_row= sel_genes,
                 column_reorder_levels = list (format_conf$branch_order, format_conf$cell_order),
                 column_legend_labels= c('branch', 'cell type'), 
                 row_legend_labels='WGCNA clusters',
                 cluster_rows=T, heat_name='norm count', center_scale=T,
                 show_column_bars = c(F, T),
                 group_order = order (module_data$MGP_PT), automatic=F)

# ----------hub score----------
hub_df <- hub_score (fil_vivo, 3, save_dir)

save_dir1 <- paste (root, 'manuscript/figure1', sep='/')
markers <- find_DE_genes (invivo, save_dir1, group.by='broad_type', label='all_vivo')
save_hub_score (markers, hub_df, 'GC10', 'EVT', save_dir)
save_hub_score (markers, hub_df, 'GC3', 'STB', save_dir)
save_hub_score (markers, hub_df, c('GC4', 'GC6'), 'CTB', save_dir)
