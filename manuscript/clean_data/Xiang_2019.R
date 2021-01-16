setwd ('utils')
library (modules)
library (Seurat)
library (tidyverse)
EM <- modules::use ('expr_mat.R')

# ----------
# Xiang 2019
# ----------

root <- 'data/Xiang_2019'
save_robj <- 'Xiang_R.Robj'

sc_data <- read.csv (paste (root, 'featurecountsLi.csv', sep='/'))
label_data <- read.csv (paste (root, 'humanIds.csv', sep='/'))
label_data2 <- read.csv (paste (root, 'Human_invitro2.csv', sep='/'))
label_data$date <- label_data2$Day [ match ( label_data$ID, label_data2$File ) ]

# =======inspect data characteristics=======
dim (sc_data)
# Convert calender month names to gene names
sc_data$Geneid<- EM$date_to_gene (as.character(sc_data$Geneid))

# ======= Clean count data =======
sc_counts <- sc_data
rownames(sc_counts) <- make.unique(sc_counts$Geneid)

sc_counts %>% dplyr::select (!Length & !Geneid) -> sc_counts

sc_counts [1:3, 1:4]

rownames(label_data) <- label_data$Filename
label_data %>% dplyr::select (!Filename) -> label_data
dataset <- CreateSeuratObject(counts = sc_counts, meta.data= label_data)
save (dataset, file = paste (root, save_robj, sep='/'))

# save dataset for GPLVM in python
markers <- read.csv (paste (root, 'marker_cluster.csv', sep='/'))
load (paste (root, save_robj, sep='/'))

top_DE <- markers %>% group_by(cluster) %>% top_n(n = 100, wt = avg_logFC)
selected_cluster <- c(0, 2, 6)
DE_genes <- top_DE [ top_DE$cluster %in% selected_cluster, ]
assay_data <- data.frame (dataset[['RNA']]@data)
selected_data <- assay_data [, dataset$seurat_clusters %in% selected_cluster]
write.csv (selected_data[rownames (selected_data) %in% DE_genes$gene, ], paste (root, 'Xiang_lognorm.csv', sep='/'))
write.csv (as.character(dataset$seurat_clusters) [dataset$seurat_clusters %in%
           selected_cluster], paste (root, 'Xiang_lognorm_cluster.csv', sep='/'))
