# This script contains a list of important marker genes for cells and/or
# transcriptional factors 

#' Load the corrected label of cells
#'
#' @param metadata metadata to be added to seurat
#' @param root_dir directory to the correct cell labels
#' @param file_name file_name of the correct cell labes
#' @param sheetName the excel sheet name
#' @param barcode_col the column storing the barcode, default is the rownames
#' @param process_ID_func function to preprocess the cell barcode 
#' @export
load_correct_label <- function (metadata, root_dir, file_name, sheetName,
                                barcode_col =NULL, process_ID_func=NULL){
        correct_label <- xlsx::read.xlsx (paste (root_dir, file_name, sep='/'),
                                    sheetName=sheetName)
        colnames (correct_label) [1] <- 'ID'

        # correct the cell barcodes
        cell_IDs <- correct_label$ID
        if (!is.null (process_ID_func)){  cell_IDs <- sapply (cell_IDs,  process_ID_func)}
        correct_label$ID <- cell_IDs

        if (is.null (barcode_col)){
                order_label <- rownames (metadata)
        }else{order_label <- metadata [, barcode_col]}

        correct_reorder <- correct_label [ match (order_label, correct_label$ID),  ]
        correct_reorder %>%
                select ('ID', 'Stage', 'Revised.lineage..this.study.', 'Ident', 'Ident_2') %>%
                magrittr::set_colnames (c( 'ID', 'Stage', 'revised', 'Ident1', 'Ident2' )) -> final_meta

        final_meta <- cbind (final_meta, metadata)
        rownames (final_meta) <- rownames (metadata)
        return (final_meta)
}

#' Standardise cell labels
#'
#' @param cellname a vector of cell names
#' @param conversion the dataframe to convert cell names
#' @return a character vector based on `cellname`
#' @export
standardise_name <- function (cellname, conversion=NULL){
        # in case a factor vector is used
        x <- as.character(cellname)
        if (is.null(conversion) ) {
                data (convert_name, package='TBdev')
                conversion <- convert_name
        }
        for (i in 1:dim(conversion)[1]){
                x [x %in% as.character (conversion$original [i]) ] <-
                        as.character (conversion$to_replace [i])
        }
        return (x)
}

fill_na_labels <- function (x, reference, fill_label){
        na_field <- is.na (x@meta.data [, fill_label])
        x@meta.data [, fill_label] <- as.character (x@meta.data [, fill_label])
        x@meta.data [na_field, fill_label] <- as.character (x@meta.data [na_field, reference])
        x@meta.data [, fill_label] <- as.factor (x@meta.data[, fill_label])
        return (x)
}

#' Append string before a numeric vector
#'
#' @description Append a string A to the beginning of another string B if
#' string B can be coerced into numeric
#'
#' @param x string B
#' @param append_str string A
append_if_numeric <- function (x, append_str='D'){
        if (is.na (as.numeric (x))){ return (x) 
        }else{return (paste (append_str, x, sep='') ) }
}

#' Clean metadata of Seurat objects
#'
#' @description The cell type labels are converted into levelled factors. 'D'
#' is appended before data labels.
#' @param AP aesthetic parameters. The key one is cell_order, which contains
#' sorting cell order
#' @param cell_type_col on which columns should cell type sorting occur
#' @param date_col which column contains date information
#' @return a Seurat object
#' @export
clean_metadata <- function (x, AP=NULL, cell_type_col= c('revised', 'Type', 
                                'assigned_cluster', 'broad_type'), date_col='date'){
        AP <- return_aes_param (AP)
        x$Type <- standardise_name (x$Type)
        if (! ('revised' %in% colnames (x@meta.data))){x$revised <- x$Type}
        na_field <- x$revised %in% c('NA', NA)
        x$revised [na_field] <- x$Type [na_field]
        for (one_col in cell_type_col){
                if (one_col %in% colnames (x@meta.data) ){
                        print (paste ('releveling', one_col))
                        x@meta.data[is.na(x@meta.data [, one_col]), one_col] <- 'unknown'
                        x@meta.data[x@meta.data [, one_col] == 'NA', one_col] <- 'unknown'
                        x@meta.data [, one_col]<- partial_relevel (x@meta.data[, one_col], AP$cell_order)
                }
        }

        # merge date and cell type label
        if (date_col %in% colnames (x@meta.data)){
                # convert days into numeric figures
                x$date <- trimws(gsub ('^D', '', as.character (x$date)))

                x$date <- sapply (x$date, function (x) {append_if_numeric (x, 'D')})
                x@meta.data %>% 
                        tidyr::unite ('Type_date', c('Type', 'date'), sep='_',  remove=F) %>%
                        tidyr::unite ('revised_date', c('revised', 'date'), sep='_',  remove=F)-> x@meta.data 
                #x$date <- paste ('D', x$date, sep='') 
                x$date <- factor ( x$date, levels= gtools::mixedsort (unique (x$date)) )
                x$Type_data <- partial_relevel (x$Type_date, AP$cell_order)
                x$revised_data <- partial_relevel (x$revised_date, AP$cell_order)
        }
        return (x)
}

#' Incorporate Trophoblast stem cells into Seurat object
#'
#' @description a specific function unlikely to be helpful
#' @export
incorporate_TSC <- function (x, analyse_by='new_cluster2', add_by='revised'){
        TSC_cells <- x$revised %in% c('hTSC_OKAE', 'hTSC_TURCO') 
        TSC_data <- x [, TSC_cells]
        TSC_data$select <- 'select'
        TSC_data@meta.data [, analyse_by] <- TSC_data@meta.data [, add_by]
        x$select <- 'not select'
        merge_x <- merge_seurat (list (x, TSC_data), assays='RNA')
        return (merge_x [, merge_x@meta.data [, analyse_by] != 'hTSC'])
}

# ----------Cells----------

select_lineage <- function (x, feature, celltypes, discard=F){
        if (discard){return (x[, ! (x@meta.data [, feature] %in% celltypes) ])
        }else{return (x[, x@meta.data [, feature] %in% celltypes])}
}

# ----------dates----------

#' Convert dates into Carnegie stages
#'
#' @export
date_to_CS <- function (x, date_col = 'date'){
        data (date2CS, package='TBdev')
        new_CS <- as.character (x@meta.data [, date_col])
        for (i in 1:nrow (date2CS)){
                new_CS [ new_CS == date2CS$date [i] ] <- as.character (date2CS$CS[i])
        }
        x@meta.data [, date_col] <- new_CS
        return (x)
}
