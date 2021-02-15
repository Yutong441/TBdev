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
#' @return a factor vector
#' @importFrom magrittr %>%
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
                data (orders, package='TBdev')
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
                ori_level [regexpr (regx, ori_level) != -1] %>% mixed_order_sort () -> level_i
                level_i %>% stringr::str_extract ('^[e,i,a,u,v,s]?') %>% 
                        factor (levels=c('e', 'i', 'a', 's', 'v', 'u', '')) %>% 
                        order() -> level_order
                level_i <- level_i [level_order]
                level_i %>% stringr::str_extract ('(early|mid|late|towards)') %>% 
                        factor (levels =c('early', 'towards', 'mid', 'late') ) %>% 
                        order () -> level_order 
                level_list [[i]] <- level_i [level_order]
        }
        do.call (c, level_list) %>% as.vector () -> final_order
        final_order <- final_order [!is.na (final_order) ]
        not_mentioned <- ori_level [! (ori_level %in% final_order) ]

        if (length (final_order) != 0 & mixed_sorting == FALSE){
                x <- factor (x, levels=unique (c ( final_order , not_mentioned )) )
        }else{ x <- factor (x, levels= mixed_order_sort (unique(as.character (x) )) )
        }
        return (x)
}

#' Slight modification to `mixedsort`
#'
#' @description `gtools::mixedsort` is a wonderful function. However, if a
#' string contains a hyphen followed by a number, the function would mistake it
#' as a negative function and sort the values accordingly. This function
#' prevents `mixedsort` from doing so.
mixed_order_sort <- function (vec){
        contain_hyphen <- sum (grepl ('-[0-9]+$', vec))
        if (contain_hyphen > 0){gtools::mixedsort (vec, decreasing=T)
        }else{gtools::mixedsort (vec)}
}

#' Merge factor vectors
#' 
#' @param list_vec a list of vectors
add_level_to_factor <- function (list_vec){
        final_vec <- unlist(lapply (list_vec, as.character ))
        final_level <- unlist(lapply (list_vec, function (x){levels (factor(x) ) }))
        return (factor (final_vec, levels=final_level))
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
return_aes_param <- function (aes_param, load_format_conf=NULL){
        if (is.null (load_format_conf)) {data (format_conf, package='TBdev')
        }else{format_conf <- load_format_conf}
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

#' Subset a Seurat object
#'
#' @description This function tries to perform subsetting in a faster way by
#' subsetting only the assay to be used and the VariableFeatures and the
#' metadata. The dimred objects are not subsetted. 
#' @param x a Seurat object
#' @importFrom Seurat VariableFeatures<-
#' @export
fast_subset <- function (x, sel_columns=NULL, sel_rows=NULL, assay='RNA', 
                         slot_data='data'){
        # subset the expression matrix
        exp_mat <- Seurat::GetAssayData (x, assay=assay, slot=slot_data)
        if (!is.null (sel_columns) & !is.null (sel_rows)){
                exp_mat <- exp_mat [sel_rows, sel_rows]
        }else if (!is.null (sel_rows)){
                exp_mat <- exp_mat [sel_rows, ]
        }else if (!is.null (sel_columns)){
                exp_mat <- exp_mat [, sel_columns]
        }

        # subset the metadata
        meta <- x@meta.data
        if (!is.null (sel_columns)){meta <- meta [sel_columns, ]}
        x_new <- Seurat::CreateSeuratObject (exp_mat, meta.data=meta, assay=assay)
        x_new <- Seurat::SetAssayData (x_new, assay=assay, slot=slot_data, new.data=exp_mat)

        # obtain variable features
        var_fea <- Seurat::VariableFeatures(x_new)
        VariableFeatures (x_new) <- var_fea[var_fea %in% rownames (exp_mat)]
        rm (var_fea, exp_mat)
        return (x_new)
}

#' @importFrom magrittr %>%
#' @export
fast_read_df <- function (save_path){
        data.table::fread (save_path) %>% data.frame () %>%
                tibble::column_to_rownames ('V1')
}
