#' Obtain label positions
#'
#' @param further_repel if TRUE, the labels would be repelled away from the
#' data points as much as possible
DimPlot_labels <- function (dat, x_coord, y_coord, color_by, further_repel=T){
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
}

#' Reimplementation of DimPlot in Seurat for better graphic controls
#' 
#' @param x a Seurat object
#' @param size_highlight a character, logical or numeric vector specifying
#' which cells to magnify in size
#' @param highlight_font how large the highlighted cells should be
#' @param AP a list for `custom_color` and `theme_TB`
#' @param further_repel if TRUE, the labels would be repelled away from the
#' data points as much as possible
#' @importFrom grDevices rainbow
#' @importFrom ggplot2 aes aes_string
#' @author Yutong Chen
#' @export
gg_DimPlot <- function (x, feature, DR='pca', dims=c(1,2), size_highlight=NULL,
                        highlight_font=4, default_color=F, further_repel=T,
                        AP=NULL,...){
        AP <- return_aes_param (AP)
        dim_red <- x@reductions[[DR]]@cell.embeddings [, dims]
        feature_names <- as.factor(x@meta.data[, feature])  
        x_axis <- colnames (dim_red)[1]
        y_axis <- colnames (dim_red)[2]

        # create a vector for highlighting size
        if (is.logical (size_highlight) ){
                size_high <- c('non-select', 'select')[as.factor(size_highlight)]
        }else if (is.null (size_highlight)) {size_high <- rep ('non-select', ncol(x) ) 
        }else if (is.character (size_highlight)) {
                size_high <- c('non-select', 'select')[as.factor (colnames (x) %in% size_highlight)]
        }else{ #if the vector is a numeric string
                size_high <- c('non-select', 'select')[as.factor(1:ncol(x) %in% size_highlight)]
        }
        dim_red %>% as.data.frame () %>%
                tibble::add_column (feature=feature_names) %>%
                tibble::add_column (size_high =size_high) -> plot_data

        plot_label <- DimPlot_labels (plot_data, x_axis, y_axis, 'feature',
                                      further_repel=further_repel)
        # change the color of the highlighted points
        color_vec <- unique (plot_data$feature)
        names (color_vec) <- color_vec
        color_vec <- rainbow ( length (color_vec) )

        high_color <- unique (plot_data$feature [plot_data$size_high == 'select'] )
        color_vec [high_color] <- rainbow (length (high_color) )

        # determine what color scheme to use
        if (!default_color){ color_vec <- custom_color (plot_data$feature, AP, color_fill=T) }
        # plotting
        plot_ob <- plot_data %>%
                ggplot2::ggplot (aes_string (x=x_axis, y=y_axis ) ) +
                ggplot2::geom_point (aes (fill = feature, size=size_high, shape=size_high), 
                            alpha=1, color=AP$point_edge_color, stroke=1) +
                ggrepel::geom_text_repel (aes (x=x_mean, y=y_mean, label=feature, fontface=2), 
                                          data=plot_label,
                                          size=AP$point_fontsize,
                                          segment.color=NA, 
                                          family=AP$font_fam)+
                ggplot2::scale_size_manual (values = c('non-select'=AP$pointsize, 
                                              'select'=AP$pointsize*1.5), guide=F) +
                # shape code: 16=filled circle, 17 = filled triangle up
                ggplot2::scale_shape_manual (values = c('non-select'=AP$normal_shape, 
                                               'select'=AP$highlight_shape), guide=F) +
                ggplot2::labs (fill= '') 
        plot_ob + theme_TB ('dim_red', plot_ob=plot_ob, feature_vec=
                            plot_data$feature, color_fill=T,
                            aes_param=AP, ...) 
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

#' Produce 3D scatter plot using gg3D
#' 
#' @param plot_data a dataframe with all the plotting information
#' @param x the column containing x coordinates, similarly for `y` and `z`
#' @param all_theta the degree of azimuthal (horizontal) rotation
#' @param all_phi the degree of vertical rotation
#' @param show_axes whehther to show 3D axes
#' @param show_arrow whether to show 3D segments (arrows)
#' @param show_label whether to show text labels
#' @param label_col which column contains the label information
#' @param num_col number of columns for `facet_wrap`
#' @param axis_length the proportion that the arrow axis occupies on the entire
#' axis
#' @param lab_just adjust the labels of the arrow axes radially. You may supply
#' one value for all 3 axes, all 3 values for x, y and z axes respectively
#' @param vert_just adjust the labels of the arrow axes vertically
#' @param further_repel if TRUE, the labels would be repelled away from the
#' data points as much as possible
#' @param force_repel extent of repulsion
#' @param AP aesthetic parameters controlling arrow appearance
#' @importFrom ggplot2 aes aes_string layer Stat
#' @author Yutong Chen
#' @export
dim_red_3D <- function (plot_data, x, y, z, color_by, all_theta=0, all_phi=0,
                        show_axes=F, show_arrow=T, show_label=T, label_col=NULL,
                        num_col=NULL, axis_length=0.2, lab_just=0.05,
                        vert_just=0., hor_just=0., further_repel=F, repel_force=1, AP=NULL){
        # deal with multiple colors
        AP <- return_aes_param (AP)
        if (length (color_by) > 1 ){
                plot_data %>% reshape2::melt (measure.vars = color_by) -> plot_data
                plot_data$variable<- partial_relevel (plot_data$variable, AP$cell_order)
                color <- 'value'
        }else{color <- color_by
        }
        ggplot2::ggplot (plot_data, aes_string (x=x, y=y, z=z) ) +
                gg3D::stat_3D (aes_string (fill=color), geom='point',
                         theta=all_theta, phi=all_phi, color= AP$point_edge_color,
                         shape=AP$normal_shape, stroke=1, size=AP$pointsize) -> plot_ob

        if (length (color_by) > 1){ plot_ob <- plot_ob + ggplot2::facet_wrap (~variable, ncol=num_col) }
        if (is.numeric (plot_data [, color]) ){
                plot_ob <- plot_ob + ggplot2::scale_fill_viridis_c () + 
                theme_dim_red (AP, color_fill=T)[[1]] 
        }else{ plot_ob + theme_TB ('no_arrow', feature_vec = plot_data [,
                                   color], color_fill=T) -> plot_ob
        }
        if (show_arrow){
                plot_ob <- plot_ob + Seg3D(theta=all_theta, phi=all_phi, common_length=axis_length) +
                        Lab3D (labs = gsub ('PT','D', c(x, y, z)), theta=all_theta, phi=all_phi,
                               common_length=axis_length+lab_just, vjust=vert_just, hjust=hor_just) 
        }else if (show_axis){
                plot_ob <- plot_ob + gg3D::axes3D(theta=all_theta, phi=all_phi) +
                        gg3D::labs_3D(labs = gsub ('PT','D', c(x, y, z)), theta=all_theta, phi=all_phi) 
        }

        if (is.null(label_col)){label_col <- color}
        if (show_label){
                print ('get text labels')
                text_scale <- get_3D_label_position (plot_data, x, y, z, label_col,
                                                     further_repel=further_repel)
                plot_ob <- plot_ob + text_3D_repel (text_scale, AP, all_theta,
                                                    all_phi, 'feature',
                                                    repel_force=repel_force)
        }
        return (plot_ob)
}

#' 3D version of DimPlot 
#' 
#' @param x a Seurat object
#' @param feature on which feature the color scheme is applied
#' @param DR which dimensionality reduction to use
#' @param ... additional features to pass to `dim_red_3D`, including
#' `show_axes`, 'all_theta' and `all_phi`
#' @author Yutong Chen
#' @export
DimPlot_3D <- function (x, feature, DR='pca', dims=c(1,2,3), ...){
        col_names <- gsub ('_', '', colnames (x@reductions[[DR]]) )
        x@reductions[[DR]]@cell.embeddings %>% as.data.frame () -> x_plot
        colnames (x_plot) <- col_names
        x_plot <- cbind (x_plot, x@meta.data)
        dim_red_3D (x_plot, col_names[dims[1]], col_names[dims[2]], col_names[dims[3]], feature, ...)
}

#' Append minimum and maximum values
#' 
#' @description Append the minimum and maximum values of `ref_data` into
#' `test_data`. This is to ensure that the same scaling to the `ref_data` will
#' be applied to `test_data` in subsequent computation
add_min_max <- function (test_data, ref_data){
        ref_min <- apply (ref_data, 2, min)
        ref_max <- apply (ref_data, 2, max)
        test_scaled <- do.call (rbind, list (test_data, ref_min, ref_max))
        return (test_scaled)
}

#' Scale data into a given range
#'
#' @param vec a numeric vector
#' @param to_range range of the output data
rescaling <- function (vec, to_range=c(0,1)){
        scaled_0_1 <- (vec - min (vec))/(max (vec) - min(vec) )
        return (scaled_0_1*(to_range[2] - to_range[1]) + to_range [1] )
}

#' @importFrom magrittr %>%
#' @noRd
get_3D_label_position <- function (test_data, tx, ty, tz, tcolor, further_repel=F){
        test_data %>% dplyr::select (c(tx, ty, tz, tcolor)) %>%
                magrittr::set_colnames (c('x_axis', 'y_axis', 'z_axis', 'feature'))  %>%
                dplyr::group_by (feature) %>%
                dplyr::summarise (x = mean(x_axis), y = mean (y_axis),
                           z = mean (z_axis)) %>% data.frame () -> plot_label

        plot_label %>% dplyr::select (c(x, y, z) ) -> test_rescaled
        add_min_max (test_rescaled, ref_data=test_data %>% 
                     dplyr::select (c(tx, ty, tz))) -> test_scaled
        label_info <- c (as.character (plot_label [, 'feature']), NA, NA)
        test_scaled$feature <- label_info

        if (further_repel){
                test_data %>% dplyr::select (dplyr::all_of (c(tx, ty, tz, tcolor))) %>%
                        magrittr::set_colnames (c('x', 'y', 'z', 'feature'))  %>%
                        dplyr::mutate (feature = rep ('', nrow (test_data) ) ) %>%
                        rbind (test_scaled) -> test_scaled
        }
        return (test_scaled)
}

#' Perform coordinate transform from 3D to 2D
#'
#' @importFrom magrittr %>%
dim_3_to_2 <- function (dat, theta, phi, axes_names=c('x', 'y', 'z')){
        pmat <- plot3D::perspbox(z=diag(2), plot=F, theta=theta, phi=phi)
        dat %>% dplyr::mutate_at (axes_names, rescaling) -> dat
        XY <- plot3D::trans3D(x = dat [,axes_names[1]],
                y = dat [,axes_names[2]], z = dat [,axes_names[3]],
                pmat = pmat) %>% data.frame()
        dat [, axes_names[1]] <- XY$x
        dat [, axes_names[2]] <- XY$y
        return (dat)
}

#' 3D version of `geom_text_repel`
#'
#' @description Same usage as `geom_text_repel` by adding it after a ggplot
#' object
#' @param dat labelling data frame
#' @param AP aesthetic parameters, only need to supply `point_fontsize`
#' @param theta angle of azimuthal rotation
#' @param phi angle of vertical rotation
#' @param axes_names the column names for the x, y, and z coordinates in `dat`
#' @param repel_force extent of repelling labels
#' @return a `geom_text_repel` layer
#' @export
text_3D_repel <- function (dat, AP, theta, phi, label_col, 
                           axes_names=c('x', 'y', 'z'), repel_force=1){
        trans_dat <- dim_3_to_2 (dat, theta, phi, axes_names)
        ggrepel::geom_text_repel (ggplot2::aes_string(x=axes_names[1],
                                                      y=axes_names[2],
                                                      label=label_col),
                                  data=trans_dat, inherit.aes=F, force=repel_force,
                                  fontface='bold', size=AP$point_fontsize) %>% list ()
}

#' Add trajectory line to 3D scatterplot
#'
#' @param plot_data the dataframe for dim red plots
#' @param px, py, pz the column names in `plot_data` that corresponds to x, y
#' and z axes
#' @param pcolor the column name in `plot_data` that provide the color for
#' points
#' @param traj_data the dataframe for the trajectory data
#' @param tx, ty, tz, tcolor similar to px, py, pz, pcolor
#' @param ... pass to `dim_red_3D`
#' @importFrom ggplot2 aes_string layer Stat
#'
#' @author Yutong Chen
#' @references 
#' \url{http://htmlpreview.github.io/?https://github.com/AckerDWM/gg3D/blob/master/gg3D-vignette.html}
#' @export
dim_red_3D_traj <- function (plot_data, px, py, pz, pcolor, traj_data, tx, ty,
                             tz, tcolor, traj_color='black', all_theta=0, all_phi=0, ...){
        # because gg3D scales everything to [0, 1]
        # To add new data on top of existing graph, it is necessary to add the
        # maximum and minimum of the existing graph to enable rescaling
        print ('rescaling axes')
        tra_scaled <- traj_data %>% tidyr::drop_na () %>% dplyr::select (c(tx, ty, tz))
        tra_ref <- plot_data %>% dplyr::select (c(px, py, pz))
        tra_scaled <- add_min_max (tra_scaled, tra_ref)

        print ('add the grouping information')
        traj_data %>% tidyr::drop_na () %>% dplyr::select (!!as.symbol (tcolor) ) -> branch_epg
        branch_epg <- c(branch_epg [, tcolor] , NA, NA)

        print ('add the color information')
        # add the color information, do not show the points labelled with NA,
        # which are only used for rescaling purpose
        lab_epg <- c(rep ('black', nrow (tra_scaled) -2), 'white', 'white')
        tra_scaled %>% as.data.frame () %>% tibble::add_column (branch = branch_epg  ) %>% 
                tibble::add_column (lab_color =  lab_epg ) -> tra_scaled

        print  ('start plotting')
        dim_red_3D (plot_data, px, py, pz, pcolor, all_theta=all_theta, all_phi=all_phi, ...) +
                gg3D::stat_3D (aes_string (group= 'branch', color= 'lab_color', x=tx, y=ty, z=tz),
                         inherit.aes=F, geom='path', theta=all_theta, phi=all_phi,
                         data=tra_scaled, size=2, linetype='dashed') +
                ggplot2::scale_color_manual (values=c('black'= traj_color, 'white'=NA))+ 
                ggplot2::guides (color=F)
}
