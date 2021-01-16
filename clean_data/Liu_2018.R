library (Seurat)

# --------
# Liu 2018
# --------

root <- 'data/Liu_2018'
save_robj <- 'Liu_R.Robj'

sc_data <- read.table (paste (root, 'GSE89497_CMP.txt', sep='/'))
label_data <- data.frame (label=colnames (sc_data))
label_data$Type2 <- sapply ( label_data$label, function(x) {strsplit (as.character (x), '_') [[1]][[2]]} )
label_data$date <- sapply ( label_data$label, function (x) {strsplit (as.character (x), 'HE')[[1]][[2]]})
label_data$date <- sapply ( label_data$date , function (x) {strsplit (as.character (x), 'W') [[1]][[1]]})
label_data$date <- 7*(as.numeric (label_data$date) - 2)
label_data$date <- sapply ( label_data$date , function (x) {paste ('D', as.character (x), sep = '') })

# remove the last digit from a string if it is numeric
# e.g. make 'CTB1' into 'CTB', whereas no changes will be made for 'CTB'
remove_last_num <- function (x){
        split_x <- strsplit (x, '')[[1]]
        if (!is.na (as.numeric (split_x [length (split_x)]))){
                split_x <- split_x [-length(split_x)]
        }
        return (paste (split_x, collapse=''))
}

label_data$Type <- sapply (label_data$Type2, remove_last_num)
rownames (label_data) <- colnames (sc_data)
dataset <- CreateSeuratObject(counts = sc_data, meta.data= label_data)
save(dataset, file = paste(root, save_robj, sep='/'))
