setwd ('utils')
library (modules)
library (Seurat)
library (tidyverse)
EM <- modules::use ('expr_mat.R')

# ---------
# Zhou 2019
# ---------

root <- 'data/Zhou_2019'
save_robj <- 'Zhou_R.Robj'
label_data <- read.csv (paste (root, 'Embryo_legend.csv', sep='/'), header=FALSE)
dim (label_data)
# [1] 5911    6
label_data [1:3,]
#                V1 V2          V3             V4 V5 V6
# 1 D8_IVC2_E2_B1_1 TE Embryo_TE_F Embryo_TE_F_D8 D8  0
# 2 D8_IVC2_E2_B1_2 TE Embryo_TE_F Embryo_TE_F_D8 D8  0
# 3 D8_IVC2_E2_B1_3 TE Embryo_TE_F Embryo_TE_F_D8 D8  0
colnames (label_data) <- c('label1', 'Type', 'label2', 'label3', 'date', 'label4')
label_data$date <- sapply (label_data$label3, function(x){strsplit(
                            as.character(x), '_')[[1]][[4]]})
rownames(label_data) <- label_data$label1

sc_data <- read.table (paste (root, 'GSE109555_All_Embryo_TPM.txt', sep='/'),
                       sep='\t', header=TRUE)
dim (sc_data)

sc_data$X <- EM$date_to_gene (as.character(sc_data$X))

rownames(sc_data) <- make.unique(sc_data$X)
sc_data %>% 
        dplyr::select (!X) -> sc_counts

sc_counts <- log1p (sc_counts/100)
dataset <- CreateSeuratObject(counts = sc_counts, meta.data= label_data)
save(dataset, file = paste(root, save_robj, sep='/'))
