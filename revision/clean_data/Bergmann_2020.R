# to use the scripts in the 'utils' folder
setwd ('utils')
library (Seurat)
library (modules)
EM <- modules::use ('expr_mat.R')

# the directory to the data
root <- 'data/Bergmann_2020'
save_robj <- 'Bergmann_R.Robj'

sc_data <- read.csv (paste (root, 'featurecounts.csv', sep='/'))
label_data <- read.csv ( paste (root, 'Key.csv', sep='/') )

dim (sc_data)
# [1] 32323  2039
dim (label_data)
# [1] 2038   35

gene_name <- sc_data$Geneid
sc_data <- sc_data [, -1]
rownames (sc_data) <- EM$date_to_gene (gene_name)
head (label_data)

table (label_data$File %in% colnames (sc_data))
length (unique (label_data$ID))

data_name <- sapply (colnames (sc_data), function (x){strsplit (x, 'Aligned.sorted')[[1]][1] })

label_data <- label_data [ match (data_name, label_data$ID) ,]
rownames (label_data) <- label_data$ID
colnames (sc_data) <- data_name

# extract cell type and label information
label_data$Annotation2 [1:10]
label_data$Type <- sapply (label_data$Annotation2, function (x){strsplit (as.character (x), '_')[[1]][1]})
label_data$Type <- sapply (label_data$Type, function (x){strsplit (as.character (x), '-')[[1]][1]})

label_data$date <- as.character (label_data$Carnegie.Stage)
table (label_data$date, label_data$Type)

# clean labels
label_data$date [is.na (label_data$date) ] <- 'no_dates'
label_data$Type <- toupper (label_data$Type)

table (label_data$Type)
ori_vec <- c('ZY', 4, 8, 'CMOR', 'ESC')
new_vec <- c('Zy', '4C', '8C', 'cMor', 'maESC')
for ( i in 1:length (ori_vec) ){
        label_data$Type [label_data$Type == ori_vec [i] ] <- new_vec [i]
}
table (label_data$Type)

dataset <- CreateSeuratObject(counts = sc_data, meta.data= label_data)
save(dataset, file = paste(root, save_robj, sep='/'))
