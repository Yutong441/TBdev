setwd ('utils')
library (modules)
library (Seurat)
library (tidyverse)
library (magrittr)
library (xlsx)
EM <- modules::use ('expr_mat.R')
ML <- modules::use ('marker_list.R')

# ----------------
# Blakeley_2015
# ----------------
# expression matrix
root <- 'data/Blakeley_2015'
save_robj <- 'Blakeley_R.Robj'
sc_counts <- read.table (paste (root, 'featurecountsNiakan.txt', sep='/'), sep='\t', header=T)
dim (sc_counts)

gene_name <- sc_counts$X
sc_counts <- sc_counts [, -1]
rownames (sc_counts) <- EM$date_to_gene (gene_name)

# label data
label_data <- read.csv (paste (root, 'Niakan_legend.csv', sep='/'), header=F)
colnames (label_data) <- c('label', 'Type')  

label_data_reorder <- label_data [match (colnames (sc_counts), label_data$label), ]
rownames (label_data_reorder) <- label_data_reorder$label

# correct label data
final_metadata <- ML$load_correct_label (label_data_reorder, root_dir, sheetName = 'Stirparo_blake_key')
final_metadata$revised <- ML$standardise_name (final_metadata$revised)
rownames (final_metadata) <- label_data_reorder$label

dataset <- CreateSeuratObject (sc_counts, meta.data = final_metadata)
save(dataset, file = paste(root, save_robj, sep='/'))

# ----------------
# Blakeley_2015 from author
# ----------------
# To test whether the results from our alignment match those from the original
# authors
root <- paste (root_dir, 'data/Blakeley_2015', sep='/')
save_robj <- 'Blakeley_R.Robj'
sc_counts <- read.table (paste (root, 'published_data.txt.gz', sep='/'), sep='\t', header=T)
gene_name <- as.character(sc_counts$Gene)

library("org.Hs.eg.db") # remember to install it if you don't have it already
gene_symbol <- mapIds(org.Hs.eg.db, keys = gene_name, keytype = "ENSEMBL", column="SYMBOL")

sc_counts <- sc_counts [, -1]
gene_symbol <- gene_symbol [!is.na (gene_symbol)]
sc_counts <- sc_counts [ match (names (gene_symbol), gene_name), ]
rownames (sc_counts) <- make.unique(gene_symbol)

label_data <- data.frame (label = colnames (sc_counts))
celltype <- gsub ('[0-9]', '', colnames (sc_counts))
label_data$Type <- gsub ('.$', '', celltype)
rownames (label_data) <- label_data$label

dataset <- CreateSeuratObject (sc_counts, meta.data = label_data)
save(dataset, file = paste(root, save_robj, sep='/'))

# ------------------------------
# Blakeley_2015 second alignment
# ------------------------------
# This is the final version
# expression matrix
root <- paste (root_dir, 'data/Blakeley_2015', sep='/')
save_robj <- 'Blakeley_R.Robj'
sc_counts <- read.table (paste (root, 'featurecountsBlakeley.txt', sep='/'), sep='\t', header=T)
dim (sc_counts)
sc_counts [1:3, 1:6]

gene_name <- sc_counts$Geneid
sc_counts <- sc_counts [, -(1:6)]
rownames (sc_counts) <- EM$date_to_gene (gene_name)
colnames (sc_counts) <- gsub ('.bam$', '', colnames(sc_counts))

# label data
label_data <- read.csv (paste (root, 'Niakan_legend.csv', sep='/'), header=F, sep=',')
colnames (label_data) <- c('cell_ID', 'Type')  
file_name <- sapply (label_data$cell_ID, function(x){strsplit (as.character (x), 'GSM') [[1]] [2]})
file_name <- sapply (file_name, function(x){strsplit (x, '_') [[1]] [1]})
label_data$label <- paste ('GSM', as.vector (file_name), sep='')

label_data_reorder <- label_data [match (colnames (sc_counts), label_data$label), ]
rownames (label_data_reorder) <- label_data_reorder$label

# correct label data
final_metadata <- ML$load_correct_label (label_data_reorder, root_dir, sheetName = 'Stirparo_blake_key', barcode_col='cell_ID')
final_metadata$revised <- ML$standardise_name (final_metadata$revised)
rownames (final_metadata) <- label_data_reorder$label

dataset <- CreateSeuratObject (sc_counts, meta.data = final_metadata)
save(dataset, file = paste(root, save_robj, sep='/'))

# -------------------------------------------------
# Correct labels from Blakeley, Yan and Petropoulos
# -------------------------------------------------
label_data <- read.xlsx (paste (root_dir,
                                'data/lineage_markers/Stirparo_full_key_new.xlsx',
                                sep='/'), sheetName='Stirparo_petro_key')

colnames (label_data)
#  [1] "Study"                                 "Embryo"                                "Cell"                                 
#  [4] "Stage"                                 "Immunosurgery..Petropolous.et.al.."    "Compacted.morula..Petropolous.et.al.."
#  [7] "Pseudotime..Petropolous.et.al.."       "Assigned.lineage..original.study."     "Revised.lineage..this.study."         
# [10] "Ident"                                 "Ident_2"                              

dim (label_data)
# Petropolous [1] 1481    9
# Blakeley [1] 28 10
# Yan: [1] 58 10

# ------------
# Naive-Primed
# ------------
root <- paste (root_dir, 'data/naive_primed', sep='/')
#save_robj <- 'Vento10X_R.Robj'

sc_counts <- read.csv (paste (root, 'counts_naive.csv', sep='/'))
dim (sc_counts)
as.character (sapply (colnames (sc_counts), function(x) { strsplit (x, '_')[[1]][1] }))
sc_counts <- read.csv (paste (root, 'counts_primed.csv', sep='/'))
dim (sc_counts)
# [1] 58394   194

sc_counts <- readRDS (paste (root, 'Human_CS7_annot_umap.rds', sep='/'))
dim (sc_counts)
# [1]  1195 57490

# ------------
# Yan 2013
# ------------
root <- paste (root_dir, 'data/Yan_2013', sep='/')
save_robj <- 'Yan_R.Robj'
sc_counts <- read.table (paste (root, 'featurecountsAllYan.txt', sep='/'), sep='\t', header=T)

sc_data <- sc_counts [, -1]
rownames(sc_data) <- EM$date_to_gene (sc_counts$Geneid)

label_data <- read.table (paste (root, 'YanIDs2.txt', sep='/'), sep='\t', header=T)
label_data <- label_data [match (colnames(sc_data), label_data$File), ]
label_data <- as.data.frame (apply (label_data, 2, as.character))
rownames (label_data) <- label_data$File
unique (label_data$Type)

date_cell <- data.frame (Type = c('Oocyte', 'Zy', '2C', '4C', '8C', 'cMor', 'Blast', 'hESC', 'hES'))
date_cell$date <- c('D0', 'D1', 'D2', 'D2', 'D3', 'D4', 'D5', 'NA', 'NA')
label_data$date <- date_cell$date [match (label_data$Type, date_cell$Type)  ]
label_data$label <- label_data$File
label_data_reorder <- label_data [ match (colnames (sc_data), label_data$label), ]

ML <- modules::use ('marker_list.R')
proc_fun <- function(x){ 
        at_gsm <- strsplit (as.character (x), 'GSM')[[1]][2] 
        gsm_str <- strsplit (at_gsm, '_')[[1]][1]
        return (paste ('GSM', gsm_str, sep=''))
}
final_metadata <- ML$load_correct_label (label_data_reorder, root_dir,
                                         sheetName = 'Stirparo_yan_key',
                                         barcode_col='ID2',
                                         process_ID_func=proc_fun)

final_metadata$revised <- ML$standardise_name (final_metadata$revised)
dataset <- CreateSeuratObject (sc_data, meta.data = final_metadata)
save(dataset, file = paste(root, save_robj, sep='/'))

# ----------------
# Petropoulos_2016
# ----------------
root <- paste (root_dir, 'data/Petropoulos_2016', sep='/')
save_robj <- 'Petropoulos_R.Robj'
sc_counts <- read.table (paste (root, 'lanner_counts.txt', sep='/'), sep='\t', header=T)
dim (sc_counts)

gene_name <- sc_counts$X
sc_counts <- sc_counts [, -1]
rownames (sc_counts) <- EM$date_to_gene (gene_name)

label_data <- read.table (paste (root, 'label_data.txt', sep='/'), sep='\t', header=T)

label_data %>% 
        select (c('Source.Name', 'Characteristics.inferred.lineage.',
                'Characteristics.developmental.stage.')) %>%
        set_colnames (c('label', 'Type', 'date')) %>%
        mutate (date = sapply (date, function(x) {
                        strsplit (as.character (x), 
                                  'embryonic day')[[1]][2]} )) -> label_data

label_data$Type <- ML$standardise_name (label_data$Type)
label_data_reorder <- label_data [match (colnames (sc_counts), label_data$label), ]
rownames (label_data_reorder) <- label_data_reorder$label

final_metadata <- ML$load_correct_label (label_data_reorder, root_dir,
                                         sheetName = 'Stirparo_petro_key')
final_metadata$revised <- ML$standardise_name (final_metadata$revised)
dataset <- CreateSeuratObject (sc_counts, meta.data = final_metadata)
save(dataset, file = paste(root, save_robj, sep='/'))
