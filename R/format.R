#' Obtain the path of dataset
#'
#' @description This function returns a list of directores with `root` meaning
#' data directory and `save_dir` result directory
#' 
#' @param root_dir where the project is located. Under this root directory,
#' there must be two other directories called 'data' and 'results'
#' @param study name of the dataset to be loaded.  
#' @export
get_path  <- function (root_dir, study){
        root <- paste (root_dir, 'data', study, sep='/')
        save_dir <- paste (root_dir, 'results', study, sep='/')
        author <- strsplit (study, '_')[[1]][1]
        save_robj <- paste (author, 'R.Robj', sep='_')
        return ( list (root=root, save_dir=save_dir, robj = save_robj) )
}

override_legend_symbol <- function (AP, color_fill=T){
        if (color_fill){
                return (list(
                      ggplot2::guides( fill= guide_legend(override.aes = list(
                                        size=AP$legend_point_size, alpha=1,
                                        shape=AP$normal_shape))
                      )))
        }else{
                return (list(
                      ggplot2::guides( color= guide_legend(override.aes = list(
                                        size=AP$legend_point_size, alpha=1))
                      )))
        }
}

#' Theme for DimPlot
#'
#' @description For dimensionality reduction plots, use a minimalistic theme
#' without any borders, grids or axis.
#' @return a list containing the theme setting that can be directory
#' concatenated with ggplot object
#'
#' @importFrom ggplot2 element_text element_blank guide_legend
#' @examples
#' ggplot (x, aes (x, y) ) +
#' geom_point () +
#' theme_dim_red ()
theme_dim_red <- function (aes_param, color_fill){
        list (ggplot2::theme(panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank(), 
              panel.background = element_blank(),
              panel.border = element_blank(),
              axis.text.x = element_blank(),
              axis.text.y = element_blank(),
              axis.title.x = element_blank(),
              axis.title.y = element_blank(),
              axis.ticks.x = element_blank(),
              axis.ticks.y = element_blank(),
              plot.title = element_text (hjust=0.5, size=aes_param$point_fontsize*3,
                                         family=aes_param$font_fam),
              legend.background= element_blank(),
              legend.key = element_blank(),
              strip.background = element_blank (),
              text=element_text (size=aes_param$fontsize, family=aes_param$font_fam),
              aspect.ratio=1,
              strip.text = element_text (size=aes_param$point_fontsize*3, family=aes_param$font_fam),
              legend.text = element_text (size=aes_param$fontsize, family=aes_param$font_fam),
              legend.title = element_text (size=aes_param$fontsize, family=aes_param$font_fam) ),
              override_legend_symbol (aes_param, color_fill)[[1]]
        )
}

#' @importFrom ggplot2 element_text element_blank guide_legend
theme_dotplot <- function (aes_param = list(fontsize=15, 
                                            point_fontsize=6, 
                                            font_fam='Arial',
                                            legend_point_size=5), 
                           rotation=90, color_fill=T){
        list (ggplot2::theme(
              # panel setting
              #panel.grid.major = element_line(color='grey92'), 
              #panel.grid.minor = element_line(color='grey92'), 
              panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank(), 
              panel.background = element_blank(),
              panel.border = element_blank(),

              # axis setting
              # hjust = 0.95, make sure the labels are assigned to the bottom
              # of the x axis
              axis.text.x = element_text(angle=rotation, family=aes_param$font_fam,
                                         hjust=0.95, size=aes_param$fontsize),
              axis.text.y = element_text(family=aes_param$font_fam, size=aes_param$fontsize),
              axis.title.x = element_text(family=aes_param$font_fam, size=aes_param$fontsize),
              axis.title.y = element_text(family=aes_param$font_fam, size=aes_param$fontsize),
              axis.ticks.x = element_blank(),
              axis.ticks.y = element_blank(),

              # legend setting
              plot.title = element_text (hjust=0.5, size=aes_param$point_fontsize*3, family=aes_param$font_fam),
              legend.background= element_blank(),
              legend.key = element_blank(),
              text=element_text (size=aes_param$fontsize, family=aes_param$font_fam),
              aspect.ratio=1,
              legend.text = element_text (size=aes_param$fontsize, family=aes_param$font_fam),
              legend.title = element_text (size=aes_param$fontsize, family=aes_param$font_fam) ,

              # facet setting
              strip.background = element_blank (),
              strip.text = element_text (size=aes_param$point_fontsize*3, family=aes_param$font_fam)
              ),

              override_legend_symbol (aes_param, color_fill)[[1]]
              )
}

custom_round <- function (vec, num_out=3, more_precision=0, quantile_val=0, round_updown=F){
        min_vec <- stats::quantile(vec, quantile_val, na.rm=T)
        max_vec <- stats::quantile(vec, 1-quantile_val, na.rm=T)
        if (max_vec < 1){
                nfig <- round (log10 (1/max_vec)) + more_precision
        }else{nfig <- 0}
        if (round_updown){
                min_vec <- ceiling (min_vec*10^(nfig))/10^(nfig)
                max_vec <- floor (max_vec*10^(nfig))/10^(nfig)
        }else{
                min_vec <- round (min_vec*10^(nfig))/10^(nfig)
                max_vec <- round (max_vec*10^(nfig))/10^(nfig)
        }
        return (seq (min_vec, max_vec, length.out=num_out))
}

#' Make the x or y axis only show the min, middle and max points
#'
#' @export
custom_tick <- function (vec, x_y='y', ...){
        breaking <- custom_round (vec, 3, ...)
        if (x_y == 'y'){return (ggplot2::scale_y_continuous (breaks = breaking) )
        }else {return (ggplot2::scale_x_continuous (breaks = breaking) )}
}

#' Extend the axis scale
#'
#' @param vec a vector as either the x or y coordinates of the data
#' @param x_y whether to extend x or y axis
#' @param extend_ratio_min how much to extend along the negative dimension of
#' the axis
#' @param extend_ratio_max how much to extend along the positive dimension of
#' the axis
#' @return `xlim` or `ylim`
custom_scale <- function (vec, x_y='y', extend_ratio_min=0, extend_ratio_max=0){
        min_vec <- min(vec, na.rm=T)
        max_vec <- max(vec, na.rm=T)
        min_vec <- min_vec - (max_vec - min_vec)*extend_ratio_min
        max_vec <- max_vec + (max_vec - min_vec)*extend_ratio_max
        if (x_y == 'y'){return (ggplot2::ylim (c(min_vec, max_vec)) )
        }else {return (ggplot2::xlim (c(min_vec, max_vec)) )}
}

#' @importFrom gtools mixedsort
custom_color <- function (feature_vec, aes_param = list (color_vec=NULL, 
                                  date_color_vec=NULL), remove_NA=T){
        # determine if the feature corresponds to cell types
        # NB: the `color_vec` refers to the color_vec defined above in this
        # script
        feature_vec <- feature_vec [!is.na (feature_vec) & feature_vec != 'NA']
        if (is.factor (feature_vec) & remove_NA){
                new_level <- levels (feature_vec) [levels (feature_vec) != 'NA']
                feature_vec <- factor (feature_vec, levels = new_level)
        }
        if (length(unique (feature_vec) ) == 0){proceed <- 'not'
        }else{
                match_celltypes <- unique (feature_vec) %in% names (aes_param$color_vec)
                if (mean (match_celltypes) == 1){proceed <- T
                }else{proceed <- F}
        }
        if (proceed == T) {new_color_vec <- aes_param$color_vec
        }else if (proceed == F){
                # determine if the feature corresponds to dates
                match_dates <- gsub ('D', '', feature_vec)
                # handle any exceptions: refer to `date_color_vec` above
                match_index <- !(match_dates %in% names (aes_param$date_color_vec ) )
                match_dates <- match_dates [match_index]
                if ( mean (is.na (as.numeric (match_dates))) == 0 ){
                        all_dates <- levels(feature_vec)
                        all_dates <- all_dates [!(all_dates %in% names (aes_param$date_color_vec) )]
                        color_vec_in <- grDevices::colorRampPalette (c('lightgray', 'black')
                                                                     )( length (all_dates) )

                        names (color_vec_in) <- all_dates 
                        new_color_vec <- c( as.character (color_vec_in), as.character (aes_param$date_color_vec) )
                        names (new_color_vec) <- mixedsort (c( as.character (names (color_vec_in)), 
                                               as.character (names (aes_param$date_color_vec) )  ))
                }else{ #when incomplete matches occur
                        print ('incomplete matches in color, filling the unknown types with rainbow color')
                        no_match <- mixedsort (unique (feature_vec) [!match_celltypes])
                        matched <- aes_param$color_vec [names (aes_param$color_vec) %in% unique (feature_vec)]
                        new_color_vec <- c(as.character (matched), rainbow (length (no_match) ))
                        names (new_color_vec) <- c(as.character (names(matched)), as.character (no_match) )
                        new_names <- partial_relevel (names (new_color_vec), aes_param$cell_order )
                        new_order <- order (new_names)
                        new_color_vec <- new_color_vec[new_order]
                }
        }
        if (proceed != 'not') {return (new_color_vec)}
}

add_custom_color_discrete <- function (feature_vec, aes_param, color_fill=F){
        cus_color <- custom_color (feature_vec, aes_param)
        if (color_fill){
                ggplot2::scale_fill_manual (values=cus_color, breaks=names (cus_color))
        }else{ggplot2::scale_color_manual (values=cus_color, breaks=names (cus_color))
        }
}

add_custom_color_continuous <- function (feature_vec, aes_param, color_fill=F){
        breaks <- custom_round (feature_vec, 2, more_precision=2, quantile_val=0.)
        if (color_fill){
                ggplot2::scale_fill_continuous (type=aes_param$palette, breaks=breaks)
        }else{ggplot2::scale_color_continuous (type=aes_param$palette, breaks=breaks)
        }
}

#' Customise the colors that correspond to particular features
#' 
#' @description This function works for either continuous or discrete scale
#' @param feature_vec a vector of features for which colors will be assigned
#' @param color_fill whether `scale_color_*` or `scale_fill_*` is used
#' @return a `scale_*` object
add_custom_color <- function (feature_vec, aes_param, color_fill=T){
        if (!is.numeric (feature_vec)){
                add_custom_color_discrete (feature_vec, aes_param, color_fill)
        }else{
                add_custom_color_continuous (feature_vec, aes_param, color_fill)
        }
}


#' Obtain arrow object
#'
#' @param AP aesthetic parameter determining arrow types, relevant ones are
#' 'arrow_angle', 'arrow_length', 'arrow_type'
#' @return a ggplot arrow object
get_arrow <- function (AP){
        ggplot2::arrow(angle=AP$arrow_angle, length = grid::unit(
                       AP$arrow_length, AP$arrow_length_unit), type= AP$arrow_type)
}

#' Create miniature arrow to replace axes
#'
#' @param length_ratio the ratio of the length of the arrow relative to the
#' span of the data
#' @param nudge_ratio how much the arrow labels should move away from the arrow
#' tip as a ratio of the length of the arrow
#' @param move_x move the arrow axies position by how much to the left,
#' negative value is to the right
#' @param aes_param setting for pointsize, font_fam, point_fontsize,
#' arrow_angle, arrow_length, arrow_length_unit, arrow_type, arrow_thickness,
#' arrow_linejoin
#' @param return a list of geoms for the arrow segment and axis labels
#' @importFrom ggplot2 geom_text 
#' @export
arrow_axis <- function (plot_ob, length_ratio=0.05, nudge_ratio=0., move_x=0,
                        move_y=0, reverse_x=F, reverse_y=F, 
                        aes_param = list (pointsize=3, font_fam='Arial',
                                          point_fontsize=6) 
                        ){
        plot_build <- ggplot2::ggplot_build (plot_ob)
        if (reverse_x){rx = -1}else{rx = 1}
        if (reverse_y){ry = -1}else{ry = 1}
        origin_x <- rx*plot_build$layout$panel_scales_x[[1]]$range$range[1]
        end_x <- rx*plot_build$layout$panel_scales_x[[1]]$range$range[2]
        origin_y <- ry*plot_build$layout$panel_scales_y[[1]]$range$range[1]
        end_y <- ry*plot_build$layout$panel_scales_y[[1]]$range$range[2]
        x <- gsub ('_', ' ', plot_build$plot$labels$x)
        y <- gsub ('_', ' ', plot_build$plot$labels$y)

        dist_x <- length_ratio*(end_x - origin_x)
        dist_y <- length_ratio*(end_y - origin_y)
        dist_xy <- max (dist_x, dist_y)

        df_arrow <- data.frame (x1=c(origin_x, origin_x), 
                                y1=c(origin_y, origin_y), 
                                x2=c(origin_x+dist_x, origin_x), 
                                y2=c(origin_y, origin_y+dist_y))

        # offset
        df_arrow$x1 <- df_arrow$x1 - dist_x*move_x
        df_arrow$x2 <- df_arrow$x2 - dist_x*move_x
        df_arrow$y1 <- df_arrow$y1 - dist_y*move_y
        df_arrow$y2 <- df_arrow$y2 - dist_y*move_y

        df_arrow$axis_labels <- c(x, y)
        df_arrow %>% 
                dplyr::mutate (xlabel = x2 + 0.2*dist_x*c(1, 0)  ) %>%
                dplyr::mutate (ylabel = y2 + 0.2*dist_y*c(0, 1)  ) -> df_arrow

        return (list(
                # type = 'closed' produces a solid triangle at the arrow end
                # linejoin = 'mitre' produces triangle with sharp edges
                ggplot2::geom_segment( aes(x = x1, y = y1, xend = x2, yend = y2), data = df_arrow, 
                            arrow = get_arrow (aes_param), size = aes_param$arrow_thickness, 
                            linejoin=aes_param$arrow_linejoin),

                # x axis
                geom_text (aes (x=xlabel, y=ylabel, label=axis_labels), data =
                           df_arrow[1,], nudge_y=0, nudge_x=rx*nudge_ratio*dist_x, 
                           size=aes_param$point_fontsize, hjust='left', vjust=0.5,
                           fontface='italic', angle=0, family=aes_param$font_fam),

                # y axis
                geom_text (aes (x=xlabel, y=ylabel, label=axis_labels), data =
                           df_arrow [2,], nudge_y=ry*nudge_ratio*dist_y, nudge_x=0, 
                           size=aes_param$point_fontsize, vjust=0.5, hjust='bottom',
                           fontface='italic', angle=90, family=aes_param$font_fam),
                # add a dot at the end of the arrows
                ggplot2::geom_point (x = rx*(origin_x - dist_x*move_x), 
                                     y = ry*(origin_y - dist_y*move_y), 
                                     size=aes_param$pointsize, shape=16)
        ))
}

#' Append this function at the end of a ggplot to provide the style in this
#' paper
#'
#' @param plot_type 'dim_red' or 'dotplot'. In 'dim_red', no grid lines will be
#' shown and arrow signs are added while in 'dotplot', grid lines will be shown
#' @param color_fill whether the color schemes will be applied to
#' `scale_fill_continuous`
#' @param whether to rotate the x axis labels or not
#' @param plot_ob if 'dim_red' is selected in `plot_type`, then it is necessary
#' to supply the ggplot object. This is to extract the axis information using
#' `ggplot_build`, which is then used to draw the axis arrows. I am attempting
#' to create a new geom layer but so far `ggproto` does not seem to allow me to
#' extract such information. 
#' NB: if automated coloring is desired, it is also essential to input plot_ob
#' @param rotation by how much degree the x axis text is rotated
#' @param length_ratio the length of the arrows for dim_plot as a percentage of
#' the plot size
#' @param aes_param a list of parameters controlling plot aesthetics
#' @param ... keyword argument to `arrow_axis`
#' @return a list of ggplot layers that can be concatenated to a ggplot object
#' @export
theme_TB <- function (plot_type='no_arrow', plot_ob=NULL, feature_vec=NULL,
                      color_fill=F, color_vec = NULL, rotation=90, 
                      aes_param = list(fontsize=15, point_fontsize=6,
                                       font_fam='Arial', pointsize=3,
                                       legend_point_size=5),
                      ...){
        aes_param <- return_aes_param (aes_param)
        if (plot_type == 'dim_red'){
                theme_list <- append ( theme_dim_red (aes_param, color_fill), arrow_axis (
                                         plot_ob, aes_param=aes_param, ...) )
        }
        if (plot_type == 'no_arrow'){
                theme_list <- theme_dim_red (aes_param, color_fill) 
        }
        if (plot_type == 'dotplot'){
                theme_list <- theme_dotplot (aes_param, rotation=rotation, color_fill=color_fill) 
        }
        if (!is.null (feature_vec)){
                if (is.numeric (feature_vec)){theme_list <- theme_list [1:length(theme_list) != 2]}
                theme_list <- append (theme_list, add_custom_color (feature_vec, 
                                                        aes_param, color_fill) )
        }
        return (theme_list)
}
