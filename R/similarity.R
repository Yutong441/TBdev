# -----------------------------------------------------
# calculate cell-cell similarity by Pearson correlation
# -----------------------------------------------------

#' Compute the difference between every single pair of cells
#'
#' @description Compute the distance matrix using a user-defined metric
#' @param x a Seurat object
#' @param method 'correlation', 'partial_corr', or 'distance'
#' @param select_cells1 a logical vector selecting the cells forming the rows
#' of the correlation matrix
#' @param select_cells2 similarly for the coluns
#' @param select_genes which genes to compare, default is variable genes
#' @return a correlation matrix with select_cells1 along the rows and
#' select_cells2 along the columns
#' @importFrom Seurat VariableFeatures
#' @importFrom magrittr %>%
#' @export
compute_all_cor <- function (x, method='correlation', assay='integrated',
                             select_cells1=NULL, select_cells2=NULL,
                             select_genes=NULL) {

        if (is.null (select_cells1)){select_cells1 <- seq (1, ncol(x))}
        if (is.null (select_cells2)){select_cells2 <- seq (1, ncol(x))}
        if (is.null (select_genes)) {select_genes <- VariableFeatures (x) }

        exp_mat <- as.matrix (x[[assay]]@data) [select_genes, ]
        if (method == 'correlation'){
                cell_cor <- stats::cor( exp_mat[, select_cells1], exp_mat [, select_cells2])
        }else if (method == 'distance'){
                exp_mat %>% t () %>% parallelDist::parDist () %>%
                        as.matrix () -> all_cor 
                cell_cor <- cell_cor [select_cells1, select_cells2]
        }else if (method == 'partial_corr'){
                cell_cor <- ppcor::pcor (exp_mat)
                cell_cor <- cell_cor$estimate [select_cells1, select_cells2]
        }else if (method == 'cosine'){
                cell_cor <- coop::cosine (exp_mat )
                cell_cor <- cell_cor[select_cells1, select_cells2]
        }
        rownames (cell_cor) <- colnames (exp_mat [, select_cells1])
        colnames (cell_cor) <- colnames (exp_mat [, select_cells2])
        return (cell_cor)
}

# -------------
# Visualisation
# -------------

#' Violin plot of cell-cell similarity
#'
#' @description Take the mean of a correlation matrix `all_cor` of the cell
#' types or features along its *column*
#' @param all_cor correlation matrix, with both colnames and rownmaes
#' @param metadata the metadata for the items in the correlation matrix, does
#' not need to be in the exact order or shape, but its rownames should match
#' the rownames and colnames of `all_cor`
#' @param feature cluster by which feature in the metadata
#' @param column_scale zero the mean across columns, make standard deviation 1
#' @param num_col number of columns in which the facets of the subplots will be
#' arranged
#' @param legend_col in how many columns the legend should be
#' @param AP aesthetic parameters, `cell_order` controls the orderring of
#' cells in the violin plot legend and other fields for plotting.
#' @importFrom magrittr %>%
#' @export
cell_violin <- function (all_cor, metadata, feature, column_scale=T,
                         row_scale=F, num_col=1, box_plot=F, AP=list (),
                         legend_col=1){
        AP <- return_aes_param (AP)
        if (length (feature) == 1){feature <- rep (feature, 2)}
        add_type <- metadata [match (colnames (all_cor), 
                                     rownames (metadata) ), feature[2] ]
        all_cor <- scaling (all_cor, row_scale, column_scale) # from 'SCHeat.R'
        all_cor %>% t () %>% as.data.frame () %>%
                tibble::add_column (cell_type =  add_type) %>%
                tidyr::gather ('all_cells', 'expr_val', -cell_type) %>%
                dplyr::group_by (all_cells, cell_type) %>%
                dplyr::summarise (mean_val = mean (expr_val)) %>%
                tidyr::spread (cell_type, mean_val) %>%
                as.data.frame ()-> plot_corr

        rownames(plot_corr) <- plot_corr$all_cells
        plot_corr$cell_type <- metadata [match (rownames (plot_corr), 
                                                rownames (metadata) ), feature [1]]
        plot_corr$cell_type <- partial_relevel (plot_corr$cell_type, AP$cell_order)
        plot_corr %>% dplyr::select (!all_cells) %>%
                tidyr::gather ('cell_type2', 'expr_val', -cell_type) -> plot_data
        ggplot2::ggplot (plot_data, ggplot2::aes (x=cell_type, y=expr_val, fill=cell_type) ) +
                ggplot2::facet_wrap (~cell_type2, ncol=num_col) +
                ggplot2::ylab ('scaled mean correlation')+
                theme_TB ('dotplot', feature_vec=plot_corr$cell_type, color_fill=T, aes_param=AP) +
                custom_tick (plot_data$expr_val, more_prec=1, x_y='y') -> plot_ob
        if (box_plot){
                plot_ob + ggplot2::geom_boxplot ()+
                        ggplot2::guides (fill=ggplot2::guide_legend (ncol=1, 
                                                override.aes=list(alpha=1))) -> plot_ob
        }else{
              plot_ob + ggplot2::geom_jitter (position=ggplot2::position_jitter(0.2), 
                                    shape=AP$normal_shape,
                                    color=AP$point_edge_color, stroke=0.3,
                                    size=AP$pointsize) +
                ggplot2::guides (fill=ggplot2::guide_legend(ncol=legend_col, override.aes=list( 
                                                           size=AP$legend_point_size, 
                                                           shape=AP$normal_shape) )) -> plot_ob
        }
}

#' Plot correlation matrix
#'
#' @param all_cor correlation matrix, with both colnames and rownmaes
#' @param metadata the metadata for the items in the correlation matrix, does
#' not need to be in the exact order or shape, but its rownames should match
#' the rownames and colnames of `all_cor`
#' @param feature cluster by which feature in the metadata
#' @export
cell_heat <- function (all_cor, metadata, features, column_scale=T,
                       row_scale=F, ...){
        all_cor <- scaling (all_cor, row_scale, column_scale) # from 'SCHeat.R'
        seurat_meta <- metadata [match (colnames (all_cor), rownames (metadata) ),]
        cor_seurat <- Seurat::CreateSeuratObject (all_cor, meta.data=seurat_meta)
        color_row <- rownames (cor_seurat)
        row_meta <- metadata [match (rownames (all_cor), rownames (metadata) ),]
        names (color_row) <- row_meta [, features[1] ]
        seurat_heat (cor_seurat, features[2], color_row, slot_data='counts',
                         show_row_names=F, ...)
}

#' Probability of fitting in river plot
#'
#' @param plot_data dataframe for the river plot
#' @param select_cells a character vector corresponding to the names of the
#' cells to plot, the column storing in names should be indicated in `meta_type`
#' @param selected_lines plot the vertical lines intersecting at the maximal
#' probability for which cell types. If NULL, all vertical lines will be shown.
#' To turn off this feature, input NA
#' @param time_ind which column stores the pseudotime information for
#' `plot_data`
#' @param band_thick thickness of the lines
#' @param band_trans transparency of the lines
#' @param branch_ind which branch to plot
#' @param sel_branch column of `plot_data` that contains the branch information
#' @param meta dataframe for plotting a color bar that serves as a reference of
#' the cell types that occur in each pseudotime
#' @param meta_time which column in `meta` that stores the pseudotime
#' informaion
#' @param meta_type which column in `meta` that stores the cell type
#' information
#' @param vjust how far down the colorbar should be from the minimum point of
#' similarity
#' @param thickness thickness of the colorbar
#' @param normalize_data whether to normalize the probability within a given
#' cell type
#' @importFrom magrittr %>%
#' @importFrom ggplot2 aes aes_string
#' @export
plot_prob_line <- function (prob_data, select_cells, selected_lines=NULL,
                            time_ind='x', band_thick=5, band_trans=0.5,
                            branch_ind='branch', sel_branch='branch1',
                            meta=NULL, meta_time='MGP_PT',
                            meta_type='broad_type', vjust=1.2, thickness=0.1,
                            normalize_data=F, remove_underscore=T, AP=NULL){
        AP <- return_aes_param (AP)
        if (normalize_data){
                prob_data %>% dplyr::mutate_at (select_cells, 
                                function(x){x/sum(x)}) -> prob_data
        }
        prob_data %>% dplyr::select (dplyr::all_of (c(select_cells, time_ind, branch_ind) ) )%>%
                dplyr::filter (!!as.symbol (branch_ind) == sel_branch) %>%
                tidyr::gather ('cell_type', 'prob', -!!as.symbol(time_ind), 
                               -!!as.symbol (branch_ind) ) -> plot_data
        plot_data %>% dplyr::group_by (cell_type) %>% dplyr::mutate ( max_prob = max (prob) ) %>% 
                dplyr::filter (prob== max_prob) %>% dplyr::ungroup() -> label_data
        if (!is.null (selected_lines)) {
                label_data %>% dplyr::filter (cell_type %in% selected_lines) -> vline_data
        }else{label_data -> vline_data}

        lim_x <- c( min (plot_data [, time_ind]), max (plot_data [, time_ind])  )
        # remove all underscores in names
        if (remove_underscore){
                plot_data$cell_type <- gsub ('_', '-', plot_data$cell_type)
                vline_data$cell_type <- gsub ('_', '-', vline_data$cell_type)
        }
        ylabel <- 'probability'
        if (normalize_data){ylabel <- paste ('normalised', ylabel)}
        ggplot2::ggplot (plot_data) +
                ggplot2::geom_line (aes_string(x=time_ind, y='prob', color='cell_type'), 
                                    size=band_thick, alpha=band_trans)+
                ggrepel::geom_text_repel (aes (x=x, y=max_prob, label=cell_type), 
                                          data=label_data, size=AP$point_fontsize, show.legend=F) +
                ggplot2::geom_vline (aes (xintercept=x, color=cell_type), 
                                     data=vline_data, linetype='dashed', show.legend=F)+
                theme_TB('dotplot', feature_vec=plot_data$cell_type, rotation=0, AP=AP)+
                ggplot2::xlab ('pseudotime') + ggplot2::ylab (ylabel) + 
                custom_tick (plot_data$prob) +
                ggplot2::xlim (lim_x) + ggplot2::theme (aspect.ratio=0.5) -> plot_ob

        if (!is.null(meta)){
                min_y <- min (plot_data [, 'prob']) - vjust
                plot_ob <- plot_ob +
                        ggplot2::geom_ribbon (aes_string (x=meta_time, fill=meta_type, 
                                                          ymin=min_y, ymax=min_y+thickness),
                                                          size=AP$pointsize, data=meta) +
                        add_custom_color (feature_vec=meta[, meta_type], AP, color_fill=T)
        }
        return (plot_ob)
}

#' Relative probability
#'
#' @description calculate probability relative to the maximum value in the same
#' category
#' @importFrom magrittr %>%
#' @export
normalize_prob <- function (x_df, features){
        x_df %>% dplyr::mutate_at (features, function (x){x/max(x)}) %>%
                magrittr::set_rownames (rownames (x_df))
}
