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
#' @param x a vector of cell names
#' @export
standardise_name <- function (cellname, conversion=NULL){
        # in case a factor vector is used
        x <- as.character(cellname)
        if (is.null(conversion) ) {data (convert_name); conversion <- convert_name
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

#' Partially relevel a vector
#'
#' @description Turn a vector into a factor vector. The order of the vector is
#' not changed, only the factor level
#'
#' @param x a vector
#' @param reorder_levels a vector specifying which celltype should be in the
#' earlier positions of of the factor level. There is no need to specify a
#' string that matches the celltype name exactly. Regular expression will be
#' used to match the string names from the beginning of the strings.
#' @param mixed_sorting if set TRUE, the `mixedsort` function from gtools package
#' will be used and will override the `reorder_levels` argument
#' @param leading_char ignore the leading character when releveling
#' @importFrom magrittr %>%
#' @importFrom gtools mixedsort
#' @examples
#' library (gtools)
#' reorder_levels <- c('Zy', '4C', '8C') # 'Zy' would be shown first in  the
#' # factor level, '4C' would be the second etc
#' celltypes <- c( '4C_D2', '8C_D2', 'Zy_D1' , '4C_D2', 'cMor_D3')
#' partial_relevel ( celltypes, reorder_levels )
#' #alternatively you may use the default ordering vector
#' @export
partial_relevel <- function (x, reorder_levels=NULL, mixed_sorting=F,
                             leading_char='', skip_char=0){
        if (is.null (reorder_levels)){
                data (orders)
                reorder_levels <- orders$cell_order
        }
        ori_level <- as.character (levels (factor (x) ))
        reorder_levels <- as.character (reorder_levels)
        level_list <- list ()
        for (i in 1:length (reorder_levels)){
                if (leading_char!= ''){
                        regx <- paste ('^', '^[e,i,a,u,v,s]?', leading_char, reorder_levels[i], sep='')
                }else{ 
                        escape_char <- paste ( rep ('.', skip_char), sep='' )
                        regx <- paste ('^' , '^[e,i,a,v,s,u]?', escape_char, reorder_levels[i], sep='') }
                # regexpr can handle regular expression but not grep
                ori_level [regexpr (regx, ori_level) != -1] %>% mixedsort () -> level_i
                level_i %>% stringr::str_extract ('^[e,i,a,u,v,s]?') %>% 
                        factor (levels=c('e', 'i', 'a', 's', 'v', 'u', '')) %>% 
                        order() -> level_order
                level_list [[i]] <- level_i [level_order]
        }
        do.call (c, level_list) %>% as.vector () -> final_order
        final_order <- final_order [!is.na (final_order) ]
        not_mentioned <- ori_level [! (ori_level %in% final_order) ]

        if (length (final_order) != 0 & mixed_sorting == FALSE){
                x <- factor (x, levels=unique (c ( final_order , not_mentioned )) )
        }else{ x <- factor (x, levels= mixedsort (unique(as.character (x) )) )
        }
        return (x)
}

#' Append string before a numeric vector
#'
#' @description Append a string A to the beginning of another string B if
#' string B can be coerced into numeric
#'
#' @param x string B
#' @param append_str string A
#' @export
append_if_numeric <- function (x, append_str='D'){
        if (is.na (as.numeric (x))){ return (x) 
        }else{return (paste (append_str, x, sep='') ) }
}

#' Clean metadata of Seurat objects
#'
#' @description The cell type labels are converted into levelled factors. 'D'
#' is appended before data labels.
#' @export
clean_metadata <- function (x){
        AP <- return_aes_param (NULL)
        x$Type <- standardise_name (x$Type)
        if (! ('revised' %in% colnames (x@meta.data))){x$revised <- x$Type}
        na_field <- x$revised %in% c('NA', NA)
        x$revised [na_field] <- x$Type [na_field]
        for (one_col in c('revised', 'Type', 'assigned_cluster', 'broad_type')){
                if (one_col %in% colnames (x@meta.data) ){
                        print (paste ('releveling', one_col))
                        x@meta.data[is.na(x@meta.data [, one_col]), one_col] <- 'unknown'
                        x@meta.data[x@meta.data [, one_col] == 'NA', one_col] <- 'unknown'
                        x@meta.data [, one_col]<- partial_relevel (x@meta.data[, one_col], AP$cell_order)
                }
        }

        # merge date and cell type label
        if ('date' %in% colnames (x@meta.data)){
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

#' @export
incorporate_TSC <- function (x, analyse_by='new_cluster2'){
        TSC_cells <- x$revised %in% c('hTSC_OKAE', 'hTSC_TURCO') 
        TSC_data <- x [, TSC_cells]
        TSC_data$select <- 'select'
        TSC_data@meta.data [, analyse_by] <- TSC_data$revised
        x$select <- 'not select'
        merge_x <- Seurat::merge (x, TSC_data)
        return (merge_x [, merge_x@meta.data [, analyse_by] != 'hTSC'])
}

# ----------Cells----------

#' @export
select_lineage <- function (x, feature, celltypes, discard=F){
        if (discard){return (x[, ! (x@meta.data [, feature] %in% celltypes) ])
        }else{return (x[, x@meta.data [, feature] %in% celltypes])}
}

# ----------dates----------

#' Convert dates into Carnegie stages
#'
#' @export
date_to_CS <- function (x, date_col = 'date'){
        data (date2CS)
        new_CS <- as.character (x@meta.data [, date_col])
        for (i in 1:nrow (date2CS)){
                new_CS [ new_CS == date2CS$date [i] ] <- as.character (date2CS$CS[i])
        }
        x@meta.data [, date_col] <- new_CS
        return (x)
}

#' Obtain aesthetic settings
#' 
#' @description If the supplied list is NULL, then the default settings will be
#' returned. Otherwise, the named fields in the supplied list will replace the
#' corresponding field in the default settings
#'
#' @param aes_param a named list of aesthetic settings
#' @examples
#' # change the fontsize and font family
#' return_aes_param (list (fontsize=12, font_fam='Arial') )
#' @export
return_aes_param <- function (aes_param){
        data (format_conf)
        if (is.null (aes_param)){
                aes_param <- format_conf
        }else{
                # if there are missing fields
                matched_names <- names (format_conf) %in% names (aes_param) 
                if ( mean (matched_names) != 1  ){
                        aes_param <- append (aes_param, format_conf [names (
                                              format_conf) [!matched_names] ] )
                }
        }
        return (aes_param)
}
