#' Obtain arrow object
#'
#' @param AP aesthetic parameter determining arrow types, relevant ones are
#' 'arrow_angle', 'arrow_length', 'arrow_type'
#' @return a ggplot arrow object
get_arrow <- function (AP){
        ggplot2::arrow(angle=AP$arrow_angle, length = grid::unit(
                       AP$arrow_length, AP$arrow_length_unit), type= AP$arrow_type)
}

find_dim_name <- function (plot_build){
        dim_name <- gsub ('_', ' ', plot_build$plot$labels$x)
        dim_name <- toupper (dim_name)
        dim_name <- trimws (gsub ('[0-9]+$', '', dim_name))
        if (dim_name == 'PC'){dim_name <- 'PCA'}
        if (dim_name == 'GP'){dim_name <- 'GPLVM'}
        return (dim_name)
}

#' @importFrom magrittr %>%
get_arrow_df <- function (plot_build, rx, ry, length_ratio, x_lab, y_lab,
                          move_x, move_y, trans_ratio=0.2, div_ratio=1){
        origin_x <- rx*plot_build$layout$panel_scales_x[[1]]$range$range[1]
        end_x <- rx*plot_build$layout$panel_scales_x[[1]]$range$range[2]
        origin_y <- ry*plot_build$layout$panel_scales_y[[1]]$range$range[1]
        end_y <- ry*plot_build$layout$panel_scales_y[[1]]$range$range[2]

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

        df_arrow$axis_labels <- c(x_lab, y_lab)
        df_arrow <- arrow_df_label (df_arrow, trans_ratio, dist_x, dist_y, 
                                    div_ratio=div_ratio)
        dim_name <- find_dim_name (plot_build)
        return (list ('df'=df_arrow, 'dx'=dist_x, 'dy'=dist_y, 'ox'=origin_x,
                      'oy'=origin_y, 'mx'=move_x, 'my'=move_y, 'rx'=rx,
                      'ry'=ry, 'dim_name'=dim_name))

}

arrow_df_label <- function (df_arrow, trans_ratio, dist_x, dist_y, div_ratio){
        if (div_ratio == 1){
                df_arrow %>% 
                        dplyr::mutate (xlabel = x2 + trans_ratio*dist_x*c(1, 0)  ) %>%
                        dplyr::mutate (ylabel = y2 + trans_ratio*dist_y*c(0, 1)  ) 
        }else{
                df_arrow %>% 
                        dplyr::mutate (xlabel = (x2 + x1)/div_ratio ) %>%
                        dplyr::mutate (ylabel = (y2 + y1)/div_ratio ) 
        }
}


#' Create ggplot layers from data for plotting arrow
#'
#' @param Arr a list generated from `get_arrow_df`
#' @importFrom ggplot2 geom_text 
arrow_to_graph <- function (Arr, aes_param, nudge_ratio){
        # type = 'closed' produces a solid triangle at the arrow end
        # linejoin = 'mitre' produces triangle with sharp edges
        layer1 <- ggplot2::geom_segment( aes(x = x1, y = y1, xend = x2, yend = y2), 
                    data = Arr$df, arrow = get_arrow (aes_param), 
                    size = aes_param$arrow_thickness, 
                    linejoin=aes_param$arrow_linejoin)

        # x axis
        layer2 <- geom_text (aes (x=xlabel, y=ylabel, label=axis_labels), data =
                   Arr$df[1,], nudge_y=0, nudge_x=Arr$rx*nudge_ratio*Arr$dx, 
                   size=aes_param$point_fontsize, hjust='left', vjust=0.5,
                   fontface='italic', angle=0, family=aes_param$font_fam)

        # y axis
        layer3 <- geom_text (aes (x=xlabel, y=ylabel, label=axis_labels), data =
                   Arr$df[2,], nudge_y=Arr$ry*nudge_ratio*Arr$dy, nudge_x=0, 
                   size=aes_param$point_fontsize, vjust=0.5, hjust='bottom',
                   fontface='italic', angle=90, family=aes_param$font_fam)

        # add a dot at the end of the arrows
        #layer4 <- ggplot2::geom_point (x = Arr$rx*(Arr$ox - Arr$dx*Arr$mx), 
        #                     y = Arr$ry*(Arr$oy - Arr$dy*Arr$my), 
        #                     size=aes_param$pointsize, shape=16)

        return (list (layer1, layer2, layer3))
}

arrow_to_graph_sim <- function (Arr, aes_param, nudge_ratio, dim_elevation=0.2,
                                nudge_dimname=0){
        # type = 'closed' produces a solid triangle at the arrow end
        # linejoin = 'mitre' produces triangle with sharp edges
        layer1 <- ggplot2::geom_segment( aes(x = x1, y = y1, xend = x2, yend = y2), 
                    data = Arr$df, arrow = get_arrow (aes_param), 
                    size = aes_param$arrow_thickness, 
                    linejoin=aes_param$arrow_linejoin)
        # x axis
        layer2 <- geom_text (aes (x=xlabel, y=ylabel, label=axis_labels), data =
                   Arr$df[1,], nudge_y=0, nudge_x=Arr$rx*nudge_ratio*Arr$dx, 
                   size=aes_param$point_fontsize, hjust='center', vjust=1,
                   fontface='italic', angle=0, family=aes_param$font_fam)

        # y axis
        layer3 <- geom_text (aes (x=xlabel, y=ylabel, label=axis_labels), data =
                   Arr$df[2,], nudge_y=Arr$ry*nudge_ratio*Arr$dy, nudge_x=0, 
                   size=aes_param$point_fontsize, vjust=0, hjust='center',
                   fontface='italic', angle=90, family=aes_param$font_fam)

        # add a dot at the end of the arrows
        #layer4 <- ggplot2::geom_point (x = Arr$rx*(Arr$ox - Arr$dx*Arr$mx), 
        #                     y = Arr$ry*(Arr$oy - Arr$dy*Arr$my), 
        #                     size=aes_param$pointsize, shape=16)

        dim_lab_y <- dim_elevation*(Arr$df[2,]$ylabel - Arr$df[1,]$ylabel) + Arr$df[1,]$ylabel
        layer5 <- geom_text (aes (x=xlabel, y=dim_lab_y, label=Arr$dim_name), data =
                   Arr$df[1,], nudge_y=0, nudge_x=Arr$rx*(nudge_ratio+nudge_dimname)*Arr$dx, 
                   size=aes_param$point_fontsize*1.3, hjust='center', vjust='bottom',
                   angle=0, family=aes_param$font_fam, fontface='bold')
        return (list (layer1, layer2, layer3, layer5))
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
#' @param axis_label_sim whether to use a simple label i.e. 'dim1', 'dim2'
#' @param nudge_dimname if `axis_label_sim` is True, then how much further to
#' the right to move the dim name label
#' @param return a list of geoms for the arrow segment and axis labels
#' @export
arrow_axis <- function (plot_ob, length_ratio=0.05, nudge_ratio=0., move_x=0,
                        move_y=0, reverse_x=F, reverse_y=F, 
                        axis_label_sim=F, nudge_dimname=0,
                        aes_param = list (pointsize=3, font_fam='Arial',
                                          point_fontsize=6) 
                        ){
        plot_build <- ggplot2::ggplot_build (plot_ob)
        if (reverse_x){rx = -1}else{rx = 1}
        if (reverse_y){ry = -1}else{ry = 1}

        # extract axis labels
        if (!axis_label_sim){
                x <- gsub ('_', ' ', plot_build$plot$labels$x)
                y <- gsub ('_', ' ', plot_build$plot$labels$y)
                div_ratio <- 1
        }else{
                x <- 'dim1'; y <- 'dim2'; div_ratio <-2
                length_ratio <- length_ratio*2
        }

        Arr <- get_arrow_df (plot_build, rx, ry, length_ratio, x, y, move_x,
                             move_y, div_ratio=div_ratio)
        if (div_ratio == 1){
                return (arrow_to_graph (Arr, aes_param, nudge_ratio) )
        }else{
                return (arrow_to_graph_sim (Arr, aes_param, nudge_ratio,
                                            nudge_dimname=nudge_dimname) )
        }
}
