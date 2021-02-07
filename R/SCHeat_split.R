#' Splitting heatmap by feature
#'
#' @export
seurat_heat_split <- function (x, color_row, group.by, split.by, main_width=NULL, split_color=NULL, ...){
        heat_list <- list ()
        if (is.null (split_color)){
                split_groups <- unique (x@meta.data [, split.by])
                split_color <- rainbow ( length (split_groups) )
                names (split_color) <- split_groups
        }else{split_groups <- names (split_color) }
        for (i in 1:length(split_groups)){
                if (i==1){left_HA <- T}else{left_HA <- F}
                if (i==1){show_column_anna <- F}else{show_column_anna <- T}
                x_select <- x [, x@meta.data [, split.by] == split_groups[i] ]
                if (!is.null (main_width)) {subwidth <- main_width*ncol (x_select)/ncol(x)
                }else{subwidth<- NULL}
                heat_list [[i]] <- seurat_heat (x_select, c(group.by, split.by), color_row, 
                                                left_HA=left_HA, provided_color= split_color[i],
                                                show_column_anna = show_column_anna, 
                                                main_width=subwidth,...)
        }
        return (heat_list)
}

rep_one_param <- function (one_param, times){
        N <- length (one_param)
        if (N < times){rep_param <- c(one_param, rep (one_param[N], times-N) )
        }else{rep_param <- one_param[1:times]}
        return (rep_param)
}

rep_param_list <- function (param_list, times){
        rep_list <- lapply (param_list, function (x){rep_one_param (x, times)})
        names (rep_list) <- names (param_list)
        return (rep_list)
}

sep_one_list <- function (all_list, index){
        N <- length (all_list)
        one_list <- lapply (as.list (1:N), function(x){all_list [[x]] [index]})
        names (one_list) <- names (all_list)
        return (one_list)
}

sep_param_list <- function (param_list){
        times <- length (param_list[[1]])
        lapply (as.list (1:times), function(x){sep_one_list (param_list, x) })
}

#' Average the expression matrix of Seurat by category
#'
#' @param x a Seurat object
#' @param group.by a vector whose first item is the feature on which averaging
#' is performed. Only the first item will be used in averaging. If a second
#' item is supplied, then a metadata field will appear containing the most
#' frequent item in the second category among the first category.
#' @param select_cells on which subset of cells averaging is performed
#' @return a Seurat object
#' @importFrom magrittr %>%
#' @export
average_by_group <- function (x, group.by, color_row, select_cells=NULL, slot_data='data',
                              assay='RNA'){
        if (is.null (select_cells)){select_cells <- 1:dim(x)[2]}
        group.by1 <- group.by [1]
        meta <- x@meta.data [select_cells, ]
        cell_index <- meta [, group.by1]
        all_types <- unique (cell_index)

        # perform averaging
        sel_exp_mat <- Seurat::GetAssayData (x, slot=slot_data, assay=assay)
        sel_exp_mat <- as.matrix (sel_exp_mat [ rownames (sel_exp_mat) %in% color_row, ])
        avg_exp_mat <- lapply (as.list (all_types), function(x){ 
                        rowMeans (sel_exp_mat [, cell_index==x] ) } )
        sel_exp_mat <- do.call (cbind, avg_exp_mat)

        # clean up, create Seurat object
        colnames (sel_exp_mat) <- all_types
        all_types %>% data.frame () %>% magrittr::set_colnames (group.by1) -> meta_data
        rownames (meta_data) <- all_types
        sel_seurat <- Seurat::CreateSeuratObject (sel_exp_mat, meta.data=meta_data)

        # add metadata if possible
        if (length (group.by) > 1 ){
                group.by2 <- group.by[2]
                table (meta [, group.by2], meta[, group.by1]) %>%
                        data.frame () %>% dplyr::group_by (Var2) %>%
                        dplyr::top_n (1, wt=Freq) -> freq_table
                matching_index <- match ( sel_seurat@meta.data [, group.by1], freq_table$Var2 )
                sel_seurat@meta.data [, group.by2] <- freq_table$Var1 [matching_index]
        }
        return (sel_seurat)
}

#' Highlight parts of a heatmap
#' 
#' @description Same as \code{seurat_heat} except that a selected part of the data
#' is shown at a larger scale and the remainder a smaller scale
#' 
#' @param x a Seurat object
#' @param select_cells a logical vector specifying which part of `x` to magnify
#' @param average average the expression across non highlighted cell types
#' @param return_sep whether the highlighted or unhighlighted heatmaps will be
#' returned as a list (TRUE) or as a ComplexHeatmapList (FALSE).
#' @param row_scale whether to perform row scaling
#' @param color_scale the scale of color gradient for heatmap. By default, the
#' `heatmap_color` field of the list supplied to `AP` should define the color
#' for min, middle and max of the dataset. Otherwise, please input a function
#' created by `circlize::colorRamp2`
#' @param center_scale whether to put white colors half way between the min and
#' max of the matrix. Otherwise, white color is placed at value 0.
#' @param show_quantile if supplied, then the color gradient would not extend
#' over the entire range of the matrix. The default 0.05 means that color
#' gradient would extend from 5% to 95% percentile. This is to avoid outliers
#' skewing the dynamic range of colors.
#' @param AP a list of aesthetic parameters. The important fields include
#' `cell_order` for ordering group labels, `gfont_fam`, `fontsize` and
#' `heatmap_color`
#'
#' @param seurat_heat_params a named list with names being the arguments to
#' `seurat_heat`. Each item must only be a vector or list with length 1 or 2.
#' If the length is 2, then the first value will be for the first
#' (unhighlighted heatmap). The second will be for the highlighted heatmap.
#' @return one or a list of ComplexHeatmaps
#' @importFrom ComplexHeatmap HeatmapAnnotation Heatmap
#' @importFrom grid gpar unit
#' @importFrom stats quantile
#' @importFrom circlize colorRamp2
#' @export
seurat_heat_highlight <- function (x, select_cells, color_row, group.by, 
                                   assay='RNA', slot_data='data',
                                   average=F, return_sep=F, 
                                   row_scale=F, column_scale=F,
                                   color_scale=NULL, center_scale=T,
                                   show_quantile=0.05, AP=NULL,
                                   seurat_heat_params = list ()){

        AP <- return_aes_param(AP)
        print ('standardise color scales across heatmap')
        if (row_scale){
                x <- x [color_row, ]
                exp_mat <- Seurat::GetAssayData (x, slot=slot_data, assay=assay)
                exp_mat <- scaling (exp_mat, row_scale, column_scale)
                x <- Seurat::SetAssayData (x, slot=slot_data, assay=assay, new.data=exp_mat)
        }
        if (is.null(color_scale)){
                plot_data <- Seurat::GetAssayData (x[color_row, select_cells], 
                                                   slot=slot_data, assay=assay)
                plot_data <- as.matrix (plot_data) 
                color_grad <- determine_color_gradient (plot_data, show_quantile, center_scale)
                color_scale <- color_grad [[1]]
                break_points <- color_grad [[2]]
                rm (plot_data)
        }else{break_points <- NULL}
        if (length (group.by) >1 ){
                split_groups <- unique (x@meta.data [, group.by[2]])
                split_color <- grDevices::rainbow ( length (split_groups), start=0.4, end=0.7)
                names (split_color) <- split_groups
        }else{split_color <- NULL}


        print ('specific heatmap parameters')
        if (length (seurat_heat_params) > 0){
                all_params <- rep_param_list (seurat_heat_params, 2)
                sep_list <- sep_param_list (all_params)
        }else{ sep_list <- list (seurat_heat_params, seurat_heat_params)
        }

        print ('common parameters')
        common_list <- list (color_row=color_row, group.by=group.by,
                             assay=assay, slot_data=slot_data,
                             provided_color=split_color,
                             break_points=break_points,
                             color_scale=color_scale, 
                             AP=AP, row_scale=F, column_scale=F
        )
        if (average){sel_seurat <- average_by_group (x, group.by, color_row,
                                !select_cells, slot_data, assay)
        }else{sel_seurat <- x[, !select_cells]
        }
        param1 <- c(list(x=sel_seurat), common_list, sep_list[[1]])
        H1 <- do.call (seurat_heat, param1)
        rm (param1)

        param2 <- c(list(x=x[, select_cells]), common_list, sep_list[[2]])
        if (length (group.by)==1 ) {H2 <- do.call (seurat_heat, param2)
        }else{H2 <- do.call (seurat_heat_split, param2)}
        rm (param2)
        if (!return_sep) {return (H1+H2)
        }else{return (list (H1, H2))}
}
