# ==========Append West 2019 to the rest of the in vivo dataset==========
root_dir <- '/mnt/c/Users/Yutong/Documents/bioinformatics/reproduction/'
root <- paste (root_dir, 'results/', sep='/')
merge_dir <- paste (root, 'XLYBPZ_Dylan_dir', sep='/')
all_data <- readRDS (paste (merge_dir, 'merged_blastoid.rds', sep='/') )
West <- readRDS (paste (root_dir, 'data/West_2019/West_R.rds', sep='/'))

setwd ('..')
devtools::load_all ()
West <- run_dim_red (West, run_diff_map=F, var_scale=T, normalize=T,
                     find_var_features=T, run_umap=F)

West$seurat_clusters <- paste ('C', West$seurat_clusters, sep='')
West$date <- gsub ('E', 'D', West$Date)
plot_dim_red (West, group.by= c('seurat_clusters'), DR='pca',
              plot_type='dim_red_sim')

data (lineage_markers)
#seurat_heat (West, color_row=lineage_markers, group.by = c('seurat_clusters', 'Date'), 
#                       slot_data='data', heat_name='norm count', center_scale=T,
#                       column_legend_labels=c('cell type', 'date'), row_scale=T)

# integrate West with the rest of the dataset
data (CT)
in_vivo <- all_data [, all_data$date != 'in_vitro' ]
West$paper <- 'West_2019'
combined <- merge_seurat ( list(in_vivo, West), 'RNA')
combined <- run_dim_red (combined, run_diff_map=F, var_scale=T, normalize=F,
                     find_var_features=T, run_umap=F)

combined$seurat_clusters <- paste ('C', combined$seurat_clusters, sep='')
plot_dim_red (combined, group.by= c('paper', 'seurat_clusters'), DR='pca', plot_type='dim_red_sim')
seurat_heat (combined, color_row=lineage_markers, group.by = c('seurat_clusters', 'date'), 
                       slot_data='data', heat_name='norm count', center_scale=T,
                       column_legend_labels=c('cell type', 'date'), row_scale=T)

table (combined$final_cluster, combined$seurat_clusters)
revision_types <- c('CTB-early', 'STB-late-1', 'STB-late-2', 'TB-late',
                    'CTB-towards-EVT-2', 'TB-mid', 'TB-early',
                    'CTB-towards-STB', 'ICM-late', 'EVT-early', 'ICM-early',
                    'EPI-1', 'STB-early', 'PE', 'EVT-late',
                    'CTB-towards-EVT-3', 'CTB-towards-EVT-1', 'cleavage',
                    'EVT-early', 'CTB-towards-EVT-1', 'STB-late-3',
                    'STB-late-2', 'EPI-2')

names (revision_types) <- paste ('C', 0:22, sep='')
combined$final_cluster <- revision_types [match (combined$seurat_clusters, names (revision_types))]
plot_dim_red (combined, group.by= c('final_cluster'), DR='pca', plot_type='dim_red_sim')

revision_broad <- gsub ('-.*$', '', revision_types)
names (revision_broad) <- paste ('C', 0:22, sep='')
combined$broad_type <- revision_broad [match (combined$seurat_clusters, names (revision_broad))]
plot_dim_red (combined, group.by= c('broad_type', 'final_cluster'), DR='pca', plot_type='dim_red_sim')

#West_part <- combined [, combined$paper == 'West_2019']
#plot_dim_red (West_part, group.by= c('final_cluster', 'old_clusters'), DR='pca', plot_type='dim_red_sim')
#seurat_heat (West_part, color_row=lineage_markers, group.by = c('final_cluster', 'date'), 
#                       slot_data='data', heat_name='norm count', center_scale=T,
#                       column_legend_labels=c('cell type', 'date'), row_scale=T)

in_vitro <- all_data [, all_data$date == 'in_vitro']
dataset <- merge_seurat ( list(combined, in_vitro), 'RNA')
Seurat::VariableFeatures (dataset) <- Seurat::VariableFeatures (combined)
saveRDS (dataset, file = paste(merge_dir, 'revision_data.rds', sep='/'))
