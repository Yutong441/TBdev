#' Obtain label positions
#'
#' @param further_repel if TRUE, the labels would be repelled away from the
#' data points as much as possible
DimPlot_labels <- function (dat, x_coord, y_coord, color_by, further_repel=T){
        if (color_by %in% colnames (dat)){
                if (!is.numeric (dat [, color_by])){
                        dat %>% dplyr::select (dplyr::all_of (c(x_coord, y_coord, color_by))) %>%
                                magrittr::set_colnames (c('x_axis', 'y_axis', 'feature'))  %>%
                                dplyr::group_by (feature) %>%
                                dplyr::summarise (x_mean = mean(x_axis), y_mean = mean (y_axis)) -> mean_labels
                        if (!further_repel){return (mean_labels)
                        }else{
                        dat %>% dplyr::select (dplyr::all_of (c(x_coord, y_coord, color_by))) %>%
                                magrittr::set_colnames (c('x_mean', 'y_mean', 'feature'))  %>%
                                dplyr::mutate (feature = rep ('', nrow (dat) ) ) %>%
                                rbind (mean_labels) 
                        }
                }else{return (NULL)}
        }else{return (NULL)}
}

get_one_feature_names <- function (name, seurat_ob, assay, slot_data){
        if (name %in% rownames(seurat_ob)){
                Seurat::FetchData (seurat_ob, name, slot=slot_data) -> feature_vec
        }else{
                feature_vec <- seurat_ob@meta.data [, name, drop=F]
        }
        return (feature_vec)
}

get_feature_names <- function (name, seurat_ob, assay, slot_data){
        if (length (name) == 1){
                feature_df <- get_one_feature_names (name, seurat_ob, assay, slot_data)
        }else{
                name_list <- list()
                for (i in 1:length (name))
                        name_list[[i]]  <- get_one_feature_names (
                                        name[i], seurat_ob, assay, slot_data)
                feature_df <- do.call (cbind, name_list)
        }
        colnames (feature_df) <- 'feature'
        meta <- seurat_ob@meta.data
        meta <- meta [, colnames (meta) != 'feature']
        return (cbind (meta, feature_df))
}

highlight_shape_size <- function (AP, highlight_ratio){
        layer1 <- ggplot2::scale_size_manual (values = c('non-select'=AP$pointsize, 
                                      'select'=AP$pointsize*highlight_ratio), guide=F) 
        # shape code: 16=filled circle, 17 = filled triangle up
        layer2 <- ggplot2::scale_shape_manual (values = c('non-select'=AP$normal_shape, 
                                       'select'=AP$highlight_shape), guide=F) 
        return (list (layer1, layer2))
}

get_size_high <- function (size_highlight, num_col, col_names=NULL){
        # create a vector for highlighting size
        if (is.logical (size_highlight) ){
                size_high <- c('non-select', 'select')[as.factor(size_highlight)]
        }else if (is.null (size_highlight)) {size_high <- rep ('non-select', num_col ) 
        }else if (is.character (size_highlight) & !is.null (col_names)) {
                size_high <- c('non-select', 'select')[as.factor (col_names %in% size_highlight)]
        }else{ #if the vector is a numeric string
                size_high <- c('non-select', 'select')[as.factor(1:num_col %in% size_highlight)]
        }
        return (size_high)
}

#' Reimplementation of DimPlot in Seurat for better graphic controls
#' 
#' @param x a Seurat object
#' @param size_highlight a character, logical or numeric vector specifying
#' which cells to magnify in size
#' @param highlight_ratio how much larger the highlighted cells should be
#' @param label_col which column contains the label information
#' @param AP a list for `custom_color` and `theme_TB`
#' @param further_repel if TRUE, the labels would be repelled away from the
#' data points as much as possible
#' @param reverse_x reverse x axis direction. This is because sometimes most of
#' the contents are on the left, which would obstruct the arrow axis.
#' @param reverse_y similarly reverse y axis direction.
#' @importFrom grDevices rainbow
#' @importFrom ggplot2 aes aes_string
#' @author Yutong Chen
gg_DimPlot <- function (x, feature, DR='pca', dims=c(1,2), size_highlight=NULL,
                        highlight_ratio=1.5, label_col=NULL, further_repel=T,
                        repel_force=1, reverse_x=F, reverse_y=F, assay='RNA',
                        slot_data = 'data', AP=NULL, seg_color=NA,
                        plot_type='dim_red',...){
        AP <- return_aes_param (AP)
        dim_red <- x@reductions[[DR]]@cell.embeddings [, dims]
        feature_names <- get_feature_names (feature, x, assay, slot_data)
        x_axis <- colnames (dim_red)[1]
        y_axis <- colnames (dim_red)[2]

        size_high <- get_size_high (size_highlight, ncol(x))
        dim_red %>% as.data.frame () %>% cbind (feature_names) %>%
                tibble::add_column (size_high =size_high) %>%
                dplyr::filter (!is.na (!!as.symbol (x_axis) ) ) -> plot_data

        if (is.null (label_col)){label_col <- 'feature'}
        plot_label <- DimPlot_labels (plot_data, x_axis, y_axis, label_col,
                                      further_repel=further_repel)

        # plotting
        ggplot2::ggplot (plot_data, aes_string (x=x_axis, y=y_axis ) ) +
                ggplot2::geom_point (aes_string (fill = 'feature', size='size_high', shape='size_high'), 
                            color=AP$point_edge_color, stroke=AP$edge_stroke) +
                highlight_shape_size (AP, highlight_ratio)+
                ggplot2::labs (fill= feature) -> plot_ob

        if (!is.null(plot_label)){
                plot_ob <- plot_ob +
                ggrepel::geom_text_repel (aes (x=x_mean, y=y_mean,
                                               label=feature, fontface=2), 
                                          data=plot_label,
                                          size=AP$point_fontsize,
                                          segment.color=seg_color, force=repel_force,
                                          family=AP$font_fam)
        }
        if (reverse_x){plot_ob <- plot_ob + ggplot2::scale_x_reverse ()}
        if (reverse_y){plot_ob <- plot_ob + ggplot2::scale_y_reverse ()}
        plot_ob + theme_TB (plot_type, plot_ob=plot_ob, feature_vec=
                            plot_data$feature, color_fill=T, aes_param=AP,
                    reverse_x=reverse_x, reverse_y=reverse_y,...) 
}

#' Return multiple subplots
gg_plot_dim_red <- function (x, by_group, DR='pca', dims=c(1,2),
                             highlight_font, size_highlight=NULL,
                             return_sep=F, ...){
        all_plots <- list ()
        for (i in 1:length (by_group)){
                all_plots [[i]] <- gg_DimPlot (x, by_group[i], DR, dims, 
                                               size_highlight, highlight_font, ...)
        }
        if (!return_sep) {
                if (length (all_plots) !=1 ){
                        return (ggpubr::ggarrange (plot_list=all_plots) )
                }else{return (all_plots [[1]] ) }
        }else{return (all_plots)}
}

