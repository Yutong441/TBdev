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
              
              # axis setting
              axis.text.x = element_blank(),
              axis.text.y = element_blank(),
              axis.title.x = element_blank(),
              axis.title.y = element_blank(),
              axis.ticks.x = element_blank(),
              axis.ticks.y = element_blank(),
              legend.background= element_blank(),
              legend.key = element_blank(),
              strip.background = element_blank (),
              text=element_text (size=aes_param$fontsize, family=aes_param$font_fam),
              plot.title = element_text (hjust=0.5, size=aes_param$fontsize*1.5,
                                         family=aes_param$font_fam, face='bold'),
              aspect.ratio=1,
              strip.text = element_text (size=aes_param$point_fontsize*3, 
                                         family=aes_param$font_fam, face='bold'),
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
              legend.background= element_blank(),
              legend.key = element_blank(),
              text=element_text (size=aes_param$fontsize, family=aes_param$font_fam),
              aspect.ratio=1,
              legend.text = element_text (size=aes_param$fontsize, family=aes_param$font_fam),
              legend.title = element_text (size=aes_param$fontsize, family=aes_param$font_fam) ,

              # facet setting
              strip.background = element_blank (),
              strip.text = element_text (size=aes_param$point_fontsize*3, 
                                         family=aes_param$font_fam, face='bold'),
              plot.title = element_text (hjust=0.5, size=aes_param$fontsize*1.5, 
                                         family=aes_param$font_fam, face='bold')
              ),
              override_legend_symbol (aes_param, color_fill)[[1]]
              )
}

custom_round <- function (vec, num_out=3, more_precision=0, quantile_val=0,
                          round_updown=F, round_downup=F, no_rounding=F,
                          lower_b=NULL, upper_b=NULL){
        min_vec <- stats::quantile(vec, quantile_val, na.rm=T)
        max_vec <- stats::quantile(vec, 1-quantile_val, na.rm=T)
        if (!is.null(lower_b)) {min_vec <- pmax (min_vec, lower_b)}
        if (!is.null(upper_b)) {max_vec <- pmin (min_vec, upper_b)}
        nfig <- round (pmin (log10 (1/max_vec), 0, na.rm=T)) + more_precision
        if (!no_rounding){
                if (round_updown){
                        min_vec <- floor (min_vec*10^(nfig))/10^(nfig)
                        max_vec <- ceiling (max_vec*10^(nfig))/10^(nfig)
                }else{
                        min_vec <- round (min_vec*10^(nfig))/10^(nfig)
                        max_vec <- round (max_vec*10^(nfig))/10^(nfig)
                }
                if (round_downup){
                        min_vec <- ceiling (min_vec*10^(nfig))/10^(nfig)
                        max_vec <- floor (max_vec*10^(nfig))/10^(nfig)
                }
        }
        return (seq (min_vec, max_vec, length.out=num_out))
}

custom_labeller <- function (vec, num_out=3, min_prec=0){
        min_vec <- min (vec, na.rm=T)
        max_vec <- max (vec, na.rm=T)
        nfig <- min_prec
        seq_vec <- seq (min_vec, max_vec, length.out=num_out)
        return_vec <- round (seq_vec, nfig)
        for (i in 1:5){
                if (return_vec[1]==return_vec[length(return_vec)]){
                        nfig <- nfig + 1
                        return_vec <- round (seq_vec, nfig)
                }else{break}
        }
        return (format (return_vec, nsmall=nfig))
}

#' Make the x or y axis only show the min, middle and max points
#'
#' @export
custom_tick <- function (vec=NULL, x_y='y', more_prec=3, min_prec=1, num_out=2,...){
        breaking <- function (x){custom_round (x, more_precision=more_prec, 
                                               num_out=num_out, no_rounding=T,...)}
        labelling <- function (x){custom_labeller (x, min_prec=min_prec, num_out=num_out)}
        if (x_y == 'y'){return (ggplot2::scale_y_continuous (breaks = breaking, labels=labelling) )
        }else {return (ggplot2::scale_x_continuous (breaks = breaking, labels=labelling) )}
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

#' Customise color
#'
#' @param feature_vec a vector of factors or characters of the names to plot
#' @param aes_param aesthetic parameters. The most important of which is
#' `color_vec`, a named vector for the colors with the names corresponding to
#' those in `feature_vec`. This needs not to be an exhaustive list. Any names
#' in `feature_vec` but not in `color_vec` will be automatically colored.
#' @param remove_NA remove NA values from `feature_vec`, i.e., do not color the
#' NA values
#' @param regexp which parts of the string in `feature_vec` to remove, in order
#' for it to match `color_vec`. This setting is useful when you try to cope
#' with varying names with the same `color_vec`.
#' @importFrom gtools mixedsort
custom_color <- function (feature_vec, aes_param = list (color_vec=NULL, 
                                  date_color_vec=NULL), remove_NA=T, regexp=NULL){
        # determine if the feature corresponds to cell types
        # NB: the `color_vec` refers to the color_vec defined above in this
        # script
        if (!is.null(regexp)){
                ori_feature <- feature_vec
                feature_vec <- gsub (regexp, '', feature_vec)
        }

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
        if (!is.null(regexp)){
                match_index <- names (new_color_vec) %in% unique(feature_vec)
                match_index2 <- match (names (new_color_vec), unique(feature_vec))
                names (new_color_vec) <- as.character (names (new_color_vec))
                names (new_color_vec)[match_index] <- as.character (unique(
                                ori_feature)[match_index2][match_index])
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

add_custom_color_continuous <- function (feature_vec, aes_param, more_prec=2, color_fill=F, ...){
        breaks <- custom_round (feature_vec, 2, more_precision=more_prec, quantile_val=0., round_updown=T)
        print (breaks)
        if (color_fill){
                ggplot2::scale_fill_continuous (type=aes_param$palette, breaks=breaks, limits=breaks)
        }else{
                print (aes_param$palette)
                ggplot2::scale_color_continuous (type=aes_param$palette, breaks=breaks, limits=breaks)
        }
}

#' Customise the colors that correspond to particular features
#' 
#' @description This function works for either continuous or discrete scale
#' @param feature_vec a vector of features for which colors will be assigned
#' @param color_fill whether `scale_color_*` or `scale_fill_*` is used
#' @return a `scale_*` object
add_custom_color <- function (feature_vec, aes_param, more_prec=0.2, color_fill=T){
        if (!is.numeric (feature_vec)){
                add_custom_color_discrete (feature_vec, aes_param, color_fill)
        }else{
                add_custom_color_continuous (feature_vec, aes_param, color_fill, more_prec=more_prec)
        }
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
                      more_prec=2,
                      aes_param = list(fontsize=15, point_fontsize=6,
                                       font_fam='Arial', pointsize=3,
                                       legend_point_size=5),
                      ...){
        aes_param <- return_aes_param (aes_param)
        if (plot_type == 'dim_red'){
                theme_list <- append ( theme_dim_red (aes_param, color_fill), arrow_axis (
                                         plot_ob, aes_param=aes_param, ...) )
        }
        if (plot_type == 'dim_red_sim'){
                theme_list <- append ( theme_dim_red (aes_param, color_fill), arrow_axis (
                                         plot_ob, aes_param=aes_param, axis_label_sim=T,...) )
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
                                                        aes_param, more_prec, color_fill) )
        }
        return (theme_list)
}
