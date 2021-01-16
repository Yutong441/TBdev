setwd ('utils')
library (modules)
library (Seurat)
library (tidyverse)
EM <- modules::use ('expr_mat.R')

# ----------
# Dylan 2020
# ----------
root <- 'data/Dylan_2020'
save_robj <- 'Dylan_R.Robj'

# ----------Process expression matrix----------#
sc_data <- read.table (paste (root, 'featurecountsDylanHuman.txt', sep='/'),
                       sep='\t', header=TRUE)
dim (sc_data)
# [1] 26364   102
colnames (sc_data)[1:10]
sc_data [1:3, 1:3]

sc_counts <- sc_data [, -(1:6)]
gene_id <- EM$date_to_gene (as.character(sc_data$Geneid))
rownames(sc_counts) <- make.unique(unlist(gene_id))

# ----------Process label_data----------
# The label data is from Dr Penfold while the cell labels from Dylan
# There appears to be a discrepancy between the two
label_data <- read.xlsx (paste (root, 'DylanData_Key.xlsx', sep='/'), sheetName='Sheet1')
cell_label <- read.csv(paste (root, 'humanseq_annotations.csv',
                                     sep='/'), head=F, stringsAsFactors = F)

colnames (cell_label) <- c('ID', 'Type')
cell_label <- cell_label [match (label_data$ID, cell_label$ID), ]

label_to_cell <- function(x){
        cells <- strsplit (x, '_')[[1]][2:3]
        return (paste (cells, collapse='_'))
}
label_data$Type <- sapply (cell_label$Type, label_to_cell)
label_data$Type [grep ('ESC', label_data$Type)] <- 'ESC'

colnames (sc_counts) <- sapply (colnames (sc_counts) , function(x) {strsplit (x, 'X')[[1]][[2]]})
label_data_reorder <- label_data [match (colnames (sc_counts), label_data$File),  ]

# ----------Filter out data with low quality----------#
select_cell <- label_data_reorder$Reads > 10000 & label_data_reorder$Efficiency > 0.5
sc_counts <- sc_counts [, select_cell]
label_data_reorder <- label_data_reorder [select_cell, ]
rownames(label_data_reorder) <- colnames(sc_counts)

dataset <- CreateSeuratObject(counts = na.omit(sc_counts), assay =
                              "RNA",min.cells = 0, min.features = 0,
                              meta.data=label_data_reorder)

save(dataset, file = paste(root, save_robj, sep='/'))
write.csv (label_data_reorder, paste (root, 'corrected_cell_label.csv', sep='/'))
