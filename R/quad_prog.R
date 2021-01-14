# ---------------------
# quadratic programming
# ---------------------

#' Calculate cell-cell similarity using quadratic programming from DeConRNASeq
#'
#' @param select_types for which cell types or cell features are similarity
#' estimates computed
#' @param compare_types to which cell types or cell features are those in
#' `select_types` compared with. If NULL, it would be every types except
#' `select_types`
#' @export
get_cell_similarity <- function (x, group.by, select_types, compare_types=NULL,
                                 slot_name='data', assay_name='RNA'){
        cell_types <- x@meta.data [, group.by [1] ]
        select_cell_index <- cell_types %in% select_types

        # set up reference
        if (is.null (compare_types)){
                all_types <- unique (cell_types)
                compare_types <- all_types [! (all_types %in% select_types) ]
        }
        cell_types_com <- x@meta.data [, group.by[2] ]
        select_cell_index2 <- cell_types_com %in% compare_types

        # get expression matrix
        exp_mat <- as.matrix ( Seurat::GetAssayData (x, slot=slot_name, assay=assay_name) )
        exp_mat_sel <- exp_mat [, select_cell_index]
        exp_mat_com <- exp_mat [, select_cell_index2]
        print (paste ('analysing', nrow(exp_mat_sel), 'genes and', ncol (exp_mat_sel), 'cells')) 

        com_cell_types <- cell_types_com [select_cell_index2]

        sign_mat <- lapply ( as.list (1:length (compare_types)), function(i){
                choose_cells <- com_cell_types == compare_types [i]
                if (length (choose_cells[choose_cells]) == 1){
                        return (exp_mat_com [, choose_cells])
                }else{ return (rowMeans (exp_mat_com [, choose_cells]) )}
        })

        print ( 'obtain the signature matrix')
        sign_mat <- do.call (cbind, sign_mat)
        colnames (sign_mat) <- compare_types
        rownames (sign_mat) <- rownames (exp_mat_com)

        print ('doing quadratic programming')
        cell_sim <- DeconRNASeq::DeconRNASeq (data.frame (exp_mat_sel), data.frame (sign_mat), fig=F)
        rownames (cell_sim [['out.all']]) <- colnames (exp_mat_sel)
        return (cell_sim[['out.all']])
}

#' Plot cell-cell similarity in a line graph
#'
#' @param x a Seurat object
#' @param cell_sim cell similarity matrix of N X M, where N is the sample and M
#' is the cell type. M must be equal to the column number of `x`
#' @param group.by which feature to color
#' @param DR dimensionality reduction assays in `x`
#' @param dims which dimension to plot along the x axis
#' @param aes_param aesthetic parameters passing to `theme_TB`
#' @importFrom ggplot2 aes 
#' @export
plot_cell_similarity <- function (x, cell_sim, group.by, DR='pca', dims=1,
                                  aes_param=NULL){
        aes_param <- return_aes_param (aes_param)
        x_select <- x[, match (rownames (cell_sim), colnames (x))]
        print (dim (x_select))
        ordering <- x_select@reductions[[DR]]@cell.embeddings[,dims]
        x_label <- colnames (x_select@reductions[[DR]]@cell.embeddings)[dims]

        cell_sim %>% data.frame () %>%
                tibble::add_column (x = ordering) %>%
                tibble::add_column (feature = x_select@meta.data [, group.by]) %>%
                tidyr::gather ('Type', 'similarity', -x, -feature) %>%
                dplyr::mutate (Type = gsub ('^X', '', Type) ) %>%
                dplyr::mutate (Type = partial_relevel (Type, aes_param$cell_order) ) -> plot_data

        plot_ob <- ggplot2::ggplot (plot_data, aes (x=x, y=similarity)) +
                        ggplot2::geom_point (aes (group=feature, color=feature) ) +
                        ggplot2::facet_wrap (~Type)+
                        ggplot2::xlab (x_label) +
                        theme_TB ('dotplot' , feature_vec=plot_data$feature,
                                  rotation=0, aes_param=aes_param) 
        return (plot_ob)

}

#' Plot cell-cell similarity in a dimensionality reduction plot
#'
#' @param x a Seurat object
#' @param cell_sim cell similarity matrix of N X M, where N is the sample and M
#' is the cell type. M must be equal to the column number of `x`
#' @param group.by which feature to color
#' @param DR dimensionality reduction assays in `x`
#' @param dims which dimensions to plot along the x and y axes
#' @param aes_param aesthetic parameters passing to `theme_TB`
#' @importFrom ggplot2 aes aes_string
#' @export
dim_plot_cell_similarity  <- function (x, cell_sim, group.by, DR='pca',
                                       dims=c(1,2), aes_param=NULL){
        x_select <- x[,match (rownames (cell_sim), colnames (x))]
        dim_red <- x_select@reductions[[DR]]@cell.embeddings [, dims]
        x_axis <- colnames (dim_red)[1]
        y_axis <- colnames (dim_red)[2]
        num_col <- integer (sqrt (ncol (cell_sim) ))
        cell_sim %>% cbind (dim_red) %>% data.frame () %>%
                tidyr::gather ('Type', 'similarity', -{{x_axis}}, -{{y_axis}}) %>%
                ggplot2::ggplot (aes_string (x=x_axis, y=y_axis)) +
                        ggplot2::geom_point (aes (color=similarity)) +
                        ggplot2::facet_wrap (~Type, ncol = num_col ) -> plot_ob
        return (plot_ob + theme_TB ('dim_red', plot_ob=plot_ob, 
                                    nudge_ratio=1.4, aes_param=aes_param)+ 
                ggplot2::scale_color_continuous (type='viridis') )
}
