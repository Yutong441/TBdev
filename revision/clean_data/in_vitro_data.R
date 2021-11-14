setwd ('..')
devtools::load_all ()
data (lineage_markers)

# ----------Sheridan 2021----------
root <- "/mnt/c/users/Yutong/Documents/bioinformatics/reproduction/"
data_dir <- paste (root, 'data/Sheridan_2021', sep='')
save_robj <- 'Sheridan_R.Robj'
sc_counts <- read.table (paste (data_dir, 'featurecountsE-MTAB-10438.txt', 
                        sep='/'), sep='\t', header=T)

# clean table
gene_name <- sc_counts$Geneid
sc_counts <- sc_counts [, -(1:6)]
rownames (sc_counts) <- gene_name

# label data
meta <- read.csv (paste (data_dir, 'Meta.csv', sep='/'), header=T)
meta$Type <- gsub ('^.*_', '', meta$X)
meta$Type [meta$Type== 'SCT'] <- 'STB_Sher'
meta$Type [meta$Type== 'EVT'] <- 'EVT_Sher'
meta$Type [meta$Type== 'TSC'] <- 'TSC_Sher'
meta$Type [meta$Type== 'TOM'] <- 'Okae_Sher'
meta <- meta [meta$Type %in% c('STB_Sher', 'EVT_Sher', 'TSC_Sher', 'Okae_Sher'),]

rownames (meta) <- meta$ID
sc_counts <- sc_counts [, match (rownames (meta), colnames (sc_counts))]
dataset <- Seurat::CreateSeuratObject (sc_counts, meta.data = meta)

dataset <- run_dim_red (dataset, run_diff_map=F, var_scale=T, normalize=T,
                     find_var_features=T, run_umap=F)

plot_dim_red (dataset, group.by= c('Type'), DR='pca', plot_type='dim_red_sim')
seurat_heat (dataset, color_row=lineage_markers, group.by = c('Type'), 
                       slot_data='data', heat_name='norm count', center_scale=T,
                       column_legend_labels='Type', row_scale=T)

saveRDS (dataset, file = paste(data_dir, save_robj, sep='/'))

# ----------Cinkornpumin ----------
root <- "/mnt/c/users/Yutong/Documents/bioinformatics/reproduction/"
data_dir <- paste (root, 'data/Cinkornpumin_2020', sep='')
save_robj <- 'Cinkornpumin_R.Robj'
sc_counts <- read.table (paste (data_dir, 'featurecountsGSM4603145.txt', 
                        sep='/'), sep='\t', header=T)

# clean table
gene_name <- sc_counts$Geneid
sc_counts <- sc_counts [, -(1:6)]
rownames (sc_counts) <- gene_name

# label data
meta <- read.csv (paste (data_dir, 'Meta.csv', sep='/'), header=T)
meta$Type <- as.character (meta$X)
meta$Type [grepl ('hTSC', meta$Type)] <- 'TSC_Cink'
meta$Type [grepl ('Primed', meta$Type)] <- 'PSC_Cink'
meta <-meta [meta$Type %in% c('TSC_Cink', 'PSC_Cink'),]

rownames (meta) <- meta$ID
sc_counts <- sc_counts [, match (rownames (meta), colnames (sc_counts))]
dataset <- Seurat::CreateSeuratObject (sc_counts, meta.data = meta)

dataset <- run_dim_red (dataset, run_diff_map=F, var_scale=T, normalize=T,
                     find_var_features=T, run_umap=F)

plot_dim_red (dataset, group.by= c('Type'), DR='pca', plot_type='dim_red_sim')
seurat_heat (dataset, color_row=lineage_markers, group.by = c('Type'), 
                       slot_data='data', heat_name='norm count', center_scale=T,
                       column_legend_labels='Type', row_scale=T)

saveRDS (dataset, file = paste(data_dir, save_robj, sep='/'))

# ----------Yanagida----------
root <- "/mnt/c/users/Yutong/Documents/bioinformatics/reproduction/"
data_dir <- paste (root, 'data/Yanagida_2021', sep='')
save_robj <- 'Yanagida_R.rds'
sc_counts <- read.table (paste (data_dir, 'featurecountsGSM5234756.txt', 
                        sep='/'), sep='\t', header=T)

# clean table
gene_name <- sc_counts$Geneid
sc_counts <- sc_counts [, -(1:6)]
rownames (sc_counts) <- gene_name

# label data
meta <- read.csv (paste (data_dir, 'Meta.csv', sep='/'), header=T)
meta$Type <- as.character (meta$X)
meta$date <- stringr::str_extract (meta$Type, 'Day[0-9]+')
meta$date <- gsub ('^Day', 'D', meta$date)
rownames (meta) <- meta$ID

# cell type labels
sc_old <- readRDS (paste (root, 'results/Yanagida_2021/Yanagida_R.rds', sep='/'))
meta$Type <- sc_old$broad_type [match (meta$X, colnames (sc_old))]

sc_counts <- sc_counts [, match (rownames (meta), colnames (sc_counts))]
dataset <- Seurat::CreateSeuratObject (sc_counts, meta.data = meta)

dataset <- run_dim_red (dataset, run_diff_map=F, var_scale=T, normalize=T,
                     find_var_features=T, run_umap=F)

plot_dim_red (dataset, group.by= c('Type'), DR='pca', plot_type='dim_red_sim')
seurat_heat (dataset, color_row=lineage_markers, group.by = c('Type'), 
                       slot_data='data', heat_name='norm count', center_scale=T,
                       column_legend_labels='Type', row_scale=T)

saveRDS (dataset, file = paste(data_dir, save_robj, sep='/'))

# ----------combine all datasets----------
data_list <- c('Sheridan_2021', 'Cinkornpumin_2020', 'Yanagida_2021')
all_data <- list ()
for (i in 1:length (data_list)){
        data_name <- paste (gsub ('[0-9]+$', 'R', data_list[i]), '.rds', sep='')
        dataset <- readRDS(paste (root, 'data', data_list[i], data_name, sep='/'))
        dataset$date <- 'in_vitro'
        all_data [[i]] <- dataset
}

all_data <- load_all_data (all_data=all_data, paper_ID=data_list)

# load combined in vivo set
root <- "/mnt/c/users/Yutong/Documents/bioinformatics/reproduction/"
merge_dir <- paste (root, 'results/XLYBPZ_Dylan_dir', sep='/')
combined <- readRDS (file = paste(merge_dir, 'revision_data.rds', sep='/'))
combined <- combined [, !combined$broad_type %in% c('bEPI-YANA', 'bTB-YANA', 'bPE-YANA')]
in_vivo <- combined [, combined$date != 'in_vitro']
in_vivo <- Seurat::FindVariableFeatures (in_vivo)

all_dataset <- merge_seurat (c(list(combined), all_data), 'RNA')
Seurat::VariableFeatures (all_dataset) <- Seurat::VariableFeatures (in_vivo)
all_dataset <- run_dim_red (all_dataset, run_diff_map=F, var_scale=T,
                         normalize=F, find_var_features=F, run_umap=F,
                         select_cells=all_dataset$date != 'in_vitro')

plot_dim_red (all_dataset, group.by= c('paper'), DR='pca', return_sep=F, 
              size_high=all_dataset$date=='in_vitro')

saveRDS (all_dataset, file = paste(merge_dir, 'revision_all.rds', sep='/'))
