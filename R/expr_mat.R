# ----------------------------
# Preprocess expression matrix
# ----------------------------

#' Convert date-like strings into their gene names
#'
#' @description Convert gene names that have been accidently converted into the
#' date format in excel back to gene names
#' @export
date_to_gene <- function (x){
        x <- as.character (x)
        month_name <- c('Jan-', 'Feb-', 'Mar-', 'Apr-', 'May-', 'Jun-', 'Jul-',
                        'Aug-', 'Sep-', 'Oct-', 'Nov-', 'Dec-')
        to_gene <- function (string){
                splitted <- strsplit (string, '-')[[1]]
                splitted[[2]] <- gsub ('0', '', splitted[[2]])
                return (paste (unlist(splitted), collapse=""))
        }

        to_change <- rep (FALSE, length(x))
        for (one_month in month_name){
                month_to_change <- startsWith (x, one_month)
                to_change <- to_change | month_to_change
        }

        x[to_change] <- sapply (x[to_change], to_gene)
        return (make.unique (as.character(x)))
}

#' Remove ribosomal genes from single cell dataset
#'
#' @param x can be S3 or S4 class
#' @export
remove_ribosomal_genes <- function (x, row_is_gene=TRUE){
        if (row_is_gene){
                ribo1 <- grep ('RPL', rownames(x))
                ribo2 <- grep ('RPS', rownames(x))
                ribo3 <- grep ('RSP', rownames(x))
                all_ribo <- c(ribo1, ribo2, ribo3)
                if (length (all_ribo) >0) {return (x[!all_ribo,])
                }else{return (x)}
                return (x[!c(ribo1, ribo2, ribo3), ])
        }else{
                ribo1 <- grep ('RPL', colnames(x))
                ribo2 <- grep ('RPS', colnames(x))
                ribo3 <- grep ('RSP', colnames(x))
                all_ribo <- c(ribo1, ribo2, ribo3)
                if (length (all_ribo) >0) {return (x[, -all_ribo])
                }else{return (x)}
        }
}

# ----------------------
# Save expression matrix
# ----------------------

extract_exp_mat <- function (x, assay=NULL, slot_data='data'){
        if (is.null(assay)){assay <- Seurat::DefaultAssay (x) }
        return (attr (x[[assay]], slot_data))
}

#' Print the dimension of a Seurat object
print_dim <- function (x){
        print (paste ('a dataset of', dim (x)[2], 'cells and', dim (x)[1], 'genes'))
}


#' Save Seurat object as csv
#' 
#' @description save the expression matrix and metadata of a Seurat object
#' @param x a seurat object
#' @param directory where the results are saved
#' @param select_cells which cells to be saved
#' @param select_genes which genes to be saved, or the number of top variably
#' expressed genes
#' @param save_result if FALSE, the expression matrix is the output. If no
#' directory is supplied, this argument is set FALSE automatically
#' @export
save_to_csv <- function (x, directory=NULL, assay=NULL, slot_data='data',
                         select_cells=NULL, select_genes=NULL, label='',
                         save_result=T){
        if (is.null (directory)){save_result <- F}
        dim_x <- dim (x)
        if (is.null(select_cells)){select_cells <- 1:dim_x[2]}
        x <- x [, select_cells]
        if (is.numeric(select_genes)){
                x <- Seurat::FindVariableFeatures (x, nfeatures=select_genes)
                select_genes <- Seurat::VariableFeatures (x)
        }

        save_x <- x [select_genes, ]
        exp_mat <- extract_exp_mat (save_x, assay, slot_data)

        print_dim (save_x)
        if (save_result){
                direct_data <- paste (directory, 'data', sep='/')
                if (!dir.exists (direct_data)){dir.create (direct_data)}
                utils::write.csv (exp_mat, paste (direct_data,
                                        paste ('merged_', label, '.csv', sep=''), sep='/'))

                utils::write.csv (save_x@meta.data, paste (direct_data, 
                                paste ('merged_meta_', label, '.csv', sep=''), sep='/' ))
        }else{
                return (as.matrix (exp_mat))
        }
}

#' Sparsify a large matrix
#' 
#' @param x a data.table object
#' @param batch_size to sparsify only a subset of the dataset at a time to
#' prevent memory overload. Set a lower batch_size for larger dataset.
sparsify <- function (x, batch_size=2000){
        batches <- list ()
        for (i in seq (1, ncol(x), batch_size)){
                print (paste ('processing batch', i%/%batch_size))
                max_select <- pmin ( i+batch_size-1, ncol (x) )
                batches [[i]] <- as.sparse (subset (x, select=i:max_select ))
        }
        return (do.call (cbind, batches))
}

#' Store 10X data
#' 
#' @description Store the expression matrix as a sparse matrix for easy loading
#' with the `read10X` function in Seurat
#'
#' @param x a data.table with the column names being the cell IDs
#' @param gene_name a list of genes for the expression matrix
#' @param gene_col a column in the expression matrix
#' @export
save_for_10X <- function (x, save_dir, batch_size=2000, gene_name=NULL, gene_col=NULL){
        if (is.character(x)){
                x <- data.table::fread (x, sep='\t', header=T)
        }
        print ('save cell IDs')
        cell_ID <- colnames (x) [! (colnames (x) %in% gene_col) ]
        utils::write.table (cell_ID, paste (save_dir,  'barcodes.tsv', sep='/'), 
                     sep='\t', row.names=F, col.names=F)

        print ('save gene names')
        if (!is.null(gene_col) == 1 ){
                gene_name <- as.matrix (data.table::subset (x, select = gene_col))
                gene_name <- sapply (gene_name, function(x){ strsplit(
                                     as.character(x), '_ENSG')[[1]][1] })
                x <- x [, (gene_col):=NULL]
        }
        gene_name <- as.data.frame (make.unique (gene_name))
        rownames (gene_name) <- gene_name [,1]
        utils:write.table (gene_name, paste (save_dir,  'genes.tsv', sep='/'), 
                     sep='\t', row.names=F, col.names=F)

        print ('save sparse matrix')
        sparse_x <- sparsify (x, batch_size=batch_size)
        print (paste ('dimension of the sparse matrix is ', dim (sparse_x)[1],
                      dim (sparse_x)[2]))
        print (paste ('dimension of the original matrix is ', dim (x)[1], dim
                      (x)[2] ))
        Matrix::writeMM (sparse_x, file = paste (save_dir, 'matrix.mtx', sep='/'))
}

save_seurat <- function(x, assay_name, file_name){
        new_x <- Seurat::CreateSeuratObject (counts=x[[assay_name]]@data, meta.data=x@meta.data)
        new_x[['RNA']]@data <- x[[assay_name]]@data
        save (new_x, file=file_name)
}
