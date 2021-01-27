# This script is essential for downstream analysis because it adds key
# information to the metadata. The following lines should always be run. The
# code blocks between '----------' can be run in separate sesssions.
setwd ('..')
devtools::load_all ()
root_dir <- '/mnt/c/Users/Yutong/Documents/bioinformatics/reproduction/'
root <- paste (root_dir, 'results/', sep='/')
merge_dir <- paste (root, 'XLYBPZ_Dylan_dir', sep='/')

# ----------standardise naming for in vitro cells----------
x <- load (paste (merge_dir, 'merged.Robj', sep='/') )
all_data <- get (x)

new_label <- as.character (all_data$revised)
old_lab <- c('hESC', 'hES', 'TSC_TU', 'TSC_OK', 'ESC')
new_lab <- c('hESC_YAN', 'hESC_YAN', 'hTSC_TURCO', 'hTSC_OKAE', 'hESC')
for (i in 1:length(old_lab)){
        new_label [new_label %in% old_lab[i] ] <- new_lab [i]
}

table (new_label)
table (all_data$revised)
all_data$revised <- new_label
all_data$publish <- 'published data'
all_data$publish [all_data$paper == 'Dylan_2020'] <- 'this study'

# ----------append assigned labels----------
all_data <- all_data [, ! (all_data$revised %in% c( ML$non_emb_lineage))]
all_data <- DIR$run_dim_red (all_data, run_diff_map=F, var_scale=T,
                             normalize=F, find_var_features=T, run_umap=F)
assignment <- read.csv ('labels/cluster_assignment.csv')
all_data$assigned_cluster <- as.character (assignment$assignment [match (all_data$seurat_clusters, assignment$cluster)])
all_data$broad_type <- assignment$broad [match (all_data$seurat_clusters, assignment$cluster)]
ML <- modules::use ('marker_list.R')
all_data <- ML$clean_metadata (all_data)
PD$gg_DimPlot (all_data, 'assigned_cluster')

save(all_data, file=paste (merge_dir, 'final_merged.Robj', sep='/') )

# ----------new method: clustering in vivo only----------
x <- load (file=paste (merge_dir, 'final_merged.Robj', sep='/') )
all_data <- get (x)
data (CT)
all_data2 <- all_data [, ! (all_data$revised %in% c( CT$non_emb_lineage, CT$in_vitro_cells))]
all_data2 <- run_dim_red (all_data2, run_diff_map=F, var_scale=T,
                             normalize=F, find_var_features=T, run_umap=F)

# determine cluster assignment
all_data2$cluster <- paste ('C', all_data2$seurat_clusters, sep='')
plot_dim_red (all_data2, by_group= c('revised', 'date', 'cluster', 'broad_type'), DR='pca', all_labels=T, return_sep=F)

data (lineage_markers)
seurat_heat (all_data2, color_row=lineage_markers, group.by = 
             c('cluster', 'revised'), column_rotation=90, slot_data='data',
             column_legend_labels=c('cell type', 'date'), row_scale=T, center_scale=T)
table (all_data2$cluster, all_data2$revised)

data (vivo_clust); assignment <- vivo_clust
all_data2$assigned_cluster <- as.character (assignment$assignment [match (all_data2$seurat_clusters, assignment$cluster)])
all_data2$broad_type <- assignment$broad [match (all_data2$seurat_clusters, assignment$cluster)]
plot_dim_red (all_data2, by_group= c('revised', 'date', 'assigned_cluster', 
              'broad_type'), DR='pca', all_labels=T, return_sep=F)

all_data$assigned_cluster <- as.character (all_data2$assigned_cluster [
                                           match (colnames (all_data), colnames (all_data2) ) ])
all_data$broad_type <- as.character (all_data2$broad_type [match (
                                           colnames (all_data), colnames (all_data2) ) ])
in_vitro_cells <- all_data$date == 'in_vitro'
all_data$assigned_cluster[in_vitro_cells] <- as.character (all_data$revised [in_vitro_cells])
all_data$broad_type [in_vitro_cells] <- as.character(all_data$revised [in_vitro_cells])
all_data <- clean_metadata (all_data)

# project the in vitro clusters
all_data <- run_dim_red (all_data, run_diff_map=F, var_scale=T,
                         normalize=F, find_var_features=F, run_umap=F,
                         select_cells=!in_vitro_cells)
plot_dim_red (all_data, group.by= c('revised', 'date', 'assigned_cluster', 
                'broad_type'), DR='pca', return_sep=F, size_high=in_vitro_cells)
save(all_data, file=paste (merge_dir, 'final_merged_tb.Robj', sep='/') )

# ----------save data for pseudotime analysis ----------
x <- load (paste (merge_dir, 'final_merged_tb.Robj', sep='/') )
all_data <- get (x)

# save from cMor
all_data2 <- all_data [, ! (all_data$revised %in% c(ML$in_vitro_cells, ML$non_emb_lineage,
                   'Oocyte', 'Zy', '2C', '4C', '8C', 'PE', 'EPI', 'PSA-EPI' ))]
all_data2 <- clean_metadata (all_data2)
plot_dim_red (all_data2, by_group= c('broad_type', 'revised'), DR='pca', all_labels=T, return_sep=F)
# 4000 top variably expressed genes
sup_save_dir <- paste (root, 'manuscript/figureS2', sep='/')
save_to_csv (all_data2, sup_save_dir, select_genes= 4000)

# save all data for reference 
exp_mat <- data.table::fread ( paste (sup_save_dir, 'data/merged_.csv', sep='/'), sep=',')
var_genes <- exp_mat$V1
save_to_csv (all_data, sup_save_dir, select_genes= var_genes, label='all')

# ----------load pseudotime results----------
x <- load (paste (merge_dir, 'final_merged_tb.Robj', sep='/') )
all_data <- get (x)
sup_save_dir <- paste (root, 'manuscript/figureS2', sep='/')
# append cell likelihood estimation
pseudo_prob <- read.csv (paste (sup_save_dir, 'result/cell_likelihood.csv', sep='/'), row.names=1)
cell_types <- c ( as.character (unique (all_data$assigned_cluster)), 'hTSC_OKAE', 'hTSC_TURCO'  )
meta <- all_data@meta.data 
new_meta <- cbind (meta, pseudo_prob [match (colnames (all_data), rownames (pseudo_prob) ) , cell_types]  )

# append GPLVM reductions
GP_PT <- read.csv (paste (sup_save_dir, 'data/PT_no_prior.csv', sep='/' ), row.names=1)
new_meta <- cbind (new_meta, GP_PT [match (rownames (new_meta), rownames(GP_PT) ), ] )

# append Epigralph branch information 
data_dir <- paste (root_dir, 'GPLVM/hier_GP_tf1/', sep='/')
branch_lab <- read.csv (paste (data_dir, 'data/STREAM_branch_labels.csv', sep='/' ), row.names=1)
pseudotime <- read.csv (paste (data_dir, 'result/infer_pt_matern.csv', sep='/' ), row.names=1)
new_meta$epil_branch <- branch_lab [match (rownames (new_meta), rownames (branch_lab) ), 1]

# the python labels are 'branch0', 'branch1' etc. I need to convert them into 
# more meaningful labels
assign_branch <- c('main', 'EVT_branch', 'STB_branch')
new_meta$epil_branch <- assign_branch [as.factor (new_meta$epil_branch) ]
new_meta$epil_branch <- factor (new_meta$epil_branch, levels= assign_branch)
# MGP_PT are pseudotimes derived from mixtures of GPLVM
new_meta$MGP_PT <- pseudotime [match (rownames (new_meta), rownames (branch_lab) ), 'pt_mean']

all_data@meta.data <- new_meta

# check cell types
table (all_data@meta.data [!is.na (all_data$MGP_PT), 'broad_type'])

# set the earliest pseudotime to be 0
all_data$ori_MGP_PT <- all_data$MGP_PT
all_data$MGP_PT <- all_data$MGP_PT - min(all_data$MGP_PT, na.rm=T)
save(all_data, file=paste (merge_dir, 'final_merged_tb.Robj', sep='/') )
