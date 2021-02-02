# --------------------------------------------------------
# scripts to reproduce the graphical abstract of the paper
# --------------------------------------------------------

#' Use sigmoid function to draw the curve
#' 
#' @param x a vector
#' @param slope the steepness of the curve. The larger the value, the smoother
#' the curve becomes
#' @param off_set shift the tuning point of the sigmoid curve. The higher the
#' value, the more leftwards the curve will shift
#' @return a vector
sigmoid <- function (x, slope, off_set, elevation=0){
        return (1/( 1 +exp (-(x- off_set)/slope) ) + elevation)
}

mid_point <- function (x){(min(x) + max(x))/2}

#' Plot sigmoid curve for one branch
#'
#' @param dat a dataframe
#' @param branch_id which branch(es) to plot
#' @param branch_col which column contains the branch information
#' @param slope the steepness of the curve. The larger the value, the smoother
#' the curve becomes
#' @param off_set shift the tuning point of the sigmoid curve. The higher the
#' value, the more leftwards the curve will shift
#' @param filter_lim the min and max of x axis
#' @param time_col which column contains pseudotime information
#' @param br_sign whether to flip the curve along the x axis, -1 indicates a
#' flip
#' @param color_col which column contains grouping information
#' @param thickness the thickness of the curve
#' @param select_cells which cells in `dat` to plot
#' @param text_offset how much the labels deviate from the top of the ribbon
#' @param elevation how much the curve is vertically translated. The direction
#' is the same as `br_sign`
#' @return a list of ggplot layers
#' @importFrom ggplot2 aes_string
#' @importFrom magrittr %>%
plot_one_branch <- function (dat, branch_id, branch_col='epil_branch', 
                             slope=1, off_set=5, 
                             filter_lim=c(-Inf, Inf),
                             time_col='MGP_PT', br_sign=1,
                             color_col='broad_type', thickness=0.01,
                             select_cells=NULL, label_group=NULL,
                             text_offset=0, AP=NULL, elevation=0, 
                             time_df=NULL, time_time_col=NULL,
                             time_color_col=NULL, time_group_col=NULL,
                             branch_name=NA, repel_force=1e-3, save_dir=NULL,
                             plot_branch_names=T){
        AP <- return_aes_param (AP)
        if (is.null(select_cells)){select_cells <- 1:nrow (dat)}
        dat <- dat [select_cells, ]
        dat [dat [, branch_col] %in% branch_id, ] %>%
                dplyr::filter (!!as.symbol (time_col) > filter_lim [1] ) %>%
                dplyr::filter (!!as.symbol (time_col) < filter_lim [2] ) -> sel_dat

        # change labels to time information if required
        if (!is.null (time_df) & !is.null (time_time_col) & 
            !is.null(time_color_col) & !is.null(time_group_col)){
                time_df_sel <- time_df [time_df [,time_group_col] %in% branch_id,]
                sel_dat [, color_col] <-append_break_labels (sel_dat [, time_col], 
                                        time_df_sel, time_time_col,
                                        time_color_col, save_dir,
                                        save_label=branch_name)
                label_group <- NULL
        }else{
                sel_dat [, color_col] <- append_break_labels (sel_dat [, time_col], 
                                        sel_dat, time_col, color_col, save_dir,
                                        save_label=branch_name)
        }

        sel_dat %>% dplyr::mutate (y_ax = sigmoid (!!as.symbol (time_col), slope, 
                                               off_set, elevation) ) %>%
                dplyr::mutate (y_ax = y_ax * br_sign) %>% 
                dplyr::mutate (y_min = y_ax - thickness) %>%
                dplyr::mutate (y_max = y_ax + thickness) -> plot_data
        if (br_sign == 1){
                plot_data %>% dplyr::mutate (y_min = pmax (y_min, 0)) -> plot_data
        }else if (br_sign == -1){
                plot_data %>% dplyr::mutate (y_max = pmin (y_max, 0)) -> plot_data
        }

        if (!is.null(label_group)){ 
                sel_dat %>% dplyr::filter (!!as.symbol (color_col) %in% label_group) -> lab_dat
        }else{lab_dat <- sel_dat}

        lab_dat %>% dplyr::group_by (!!as.symbol (color_col) ) %>%
                dplyr::summarise (mean_pt = mid_point (!!as.symbol (time_col) )) %>%
                dplyr::mutate (y_max = sigmoid (mean_pt, slope, off_set )+ 
                        thickness + text_offset) %>%
                dplyr::mutate (y_max = y_max*br_sign) -> label_dat

        layer1 <- ggplot2::geom_ribbon (aes_string (x=time_col, ymin='y_min', ymax='y_max', 
                                            fill=color_col), data=plot_data, color='black') 
        layer2 <- ggrepel::geom_text_repel ( aes_string (x='mean_pt', y='y_max', label=color_col), 
                             data = label_dat, inherit.aes=F, force=repel_force, 
                             size=AP$point_fontsize, family=AP$font_fam, fontface='bold')
        plot_layers <- list (layer1, layer2)
        if (!is.na (branch_name) & plot_branch_names){
                max_t <- max(plot_data [, time_col])
                max_i <- plot_data [, time_col]==max_t
                mid_y <- (plot_data$y_min [max_i] + plot_data$y_max [max_i])/2
                layer3 <- ggplot2::annotate ('text', x=max_t, y=mid_y, label=branch_name, 
                                             hjust='right', size=AP$point_fontsize, family=AP$font_fam) 
                plot_layers <- c(plot_layers, list (layer3))
        }
        return (plot_layers)
}

#' Manually draw an x-axis
#'
#' @description Sometimes you may want to put the x axis in certain locations
#' other than the periphery of the plot as in the default ggplot option. Here I
#' only implemented the version for x axis but it would be trivial for y axis
#' as well. The underlying rationale is simple. Draw the axis and put on the
#' labels.
#' @param min_x where to draw the x axis from
#' @param max_x where the axis ends
#' @param x_name name for the x axis
#' @param AP aesthetic paameters
#' @param y_level which y coordinate
#' @param shorten_ratio shorten the x axis from `max_x`. This is because
#' sometimes the axis tip is cut off from the edge of the plot.
#' @param ticks_df a dataframe for plotting the axis ticks
#' @return a list of ggplot layers, must contain a column called `x` for the
#' locations of the ticks and the labels to put on
#' @importFrom ggplot2 aes
get_time_axis <- function (min_x, max_x, x_name, AP, y_level=0,
                           shorten_ratio=0.85, ticks_df=NULL, 
                           tick_height=5e-3){
        # draw the axis itself
        layer1 <- ggplot2::geom_segment (aes (x=min_x, y=y_level, 
                 xend=max_x*shorten_ratio, yend=0), arrow=get_arrow (AP),
                 size=AP$arrow_thickness, linejoin = AP$arrow_linejoin)

        # put the label for x axis
        layer2 <- ggplot2::annotate('text', x=max_x, 
                                    y=y_level, label=x_name, size=AP$point_fontsize,
                                    family=AP$font_fam, hjust='right') 
        axis_layer <- list(layer1, layer2)
        if (!is.null (ticks_df)){
                layer3 <- ggrepel::geom_text_repel (aes (x=x, y=y_level-tick_height, 
                                label=x), data=ticks_df, vjust='bottom')
                layer4 <- ggplot2::geom_segment (aes (x=x, xend=x, 
                                y=y_level, yend=y_level-tick_height), data=ticks_df)
                axis_layer <- c(axis_layer, list(layer3, layer4))
        }
        return (axis_layer)
}

#' Draw 2 sigmoid curves
#' 
#' @param branch1 which samples in `dat` make up the first branch
#' @param select_cells1 which cells in the first branch
#' @param time_axis whether to plot the x axis
#' @param axis_label what to label the `time_axis`
#' @param vert_trans the extent of vertical translation
#' @return a list of ggplot layers
#' @seealso `plot_one_branch`
#' @examples
#' ggplot () +
#' plot_two_branch (meta, 
#'                  select_cells1 = meta$broad_type != 'STB',
#'                  select_cells2 = meta$broad_type != 'EVT',
#'                  text_offset=0.007, off_set= 5
#' ) + ggplot2::theme (legend.position='none')
#' @importFrom ggplot2 aes
#' @export
plot_two_branch <- function (dat, branch_col = 'epil_branch', 
                             color_col = 'broad_type',
                             time_col = 'MGP_PT',
                             branch1= c('TB_stem', 'STB_branch'), 
                             branch2 = c('TB_stem', 'EVT_branch'), 
                             select_cells1 = NULL,
                             select_cells2 = NULL, AP=NULL,
                             time_axis=T, axis_label='time',
                             vert_trans=0, branch_names=c(NA,NA,NA), 
                             plot_branch_names=T,
                             thickness=0.01,...){
        AP <- return_aes_param (AP)
        r1 <- range (dat [dat [, branch_col] %in% branch1, time_col])
        r2 <- range (dat [dat [, branch_col] %in% branch2, time_col])
        min_b <- pmax (r1[1], r2[1])
        max_b <- pmin (r1[2], r2[2])
        filter_lim <- c(min_b, max_b)

        # determine redundant labels
        if (is.null(select_cells1)){select_cells1 <- 1:nrow (dat)}
        if (is.null(select_cells2)){select_cells2 <- 1:nrow (dat)}
        dat1 <- dat [select_cells1, ]
        dat2 <- dat [select_cells2, ]

        uniq1 <- unique (dat1 [ dat1 [, branch_col] %in% branch1, color_col ])
        uniq2 <- unique (dat2 [ dat2 [, branch_col] %in% branch2, color_col ])
        label_group <- uniq2 [!uniq2 %in% uniq1]

        plot_layers <- c(
                plot_one_branch (dat, branch2, time_col = time_col, 
                                 color_col = color_col, select_cells= select_cells2,
                                 branch_col=branch_col, br_sign=-1, 
                                 filter_lim=filter_lim, label_group=label_group, 
                                 AP=AP, elevation=vert_trans, 
                                 branch_name=branch_names[1], thickness=thickness,
                                 plot_branch_names=plot_branch_names,...),
                plot_one_branch (dat, branch1, time_col = time_col, 
                                 color_col = color_col, select_cells= select_cells1,
                                 branch_col=branch_col, filter_lim=filter_lim, 
                                 AP=AP, elevation=vert_trans,
                                 branch_name=branch_names[2], thickness=thickness,
                                 plot_branch_names=plot_branch_names,...),
                theme_TB ('no_arrow', feature_vec=dat[, color_col], color_fill=T, AP=AP) 
        )
        if (time_axis){
                ticks_df <- data.frame (x = 0:floor(max_b*0.9) )
                axis_layer <- get_time_axis (min_b, max_b, axis_label, AP, 
                                             ticks_df=ticks_df, shorten_ratio=0.9, 
                                             tick_height=max(dat[, time_col])*5e-4)
                plot_layers <- c(plot_layers, axis_layer) 
        }
        if (length (branch_names) >=3 & plot_branch_names){
                if (!is.na (branch_names[3]) ){
                        layer3 <- ggplot2::annotate ('text', x=min_b, y=-thickness, 
                                 size=AP$point_fontsize, family=AP$font_fam, 
                                 hjust='left', vjust='top', label=branch_names[3])
                        plot_layers <- c(plot_layers, list(layer3))
                }
        }
        return (plot_layers)
}

#' Obtain break points along a time line
#'
#' @description This function finds the interval of different periods along a
#' time line. This is not trivial, because sometimes different periods overlap
#' with each other. For aesthetic reasons however, there is significant demand
#' to show the periods as being mutually exclusive (which I strongly dislike).
#' As a result, this function will first find the min and max of time periods
#' in each group. The min is taken as the starting point. If the ending point
#' is the min of the next time period. If the min of the next time period
#' precedes that of its previous one, the starting point of the next time
#' period is taken to be the max of the previous time period.
#' @param x a dataframe of all items
#' @param time_col which column in `x` has the time information
#' @param color_col which column in `x` has the grouping information
#' @param time_order order of the grouping labels in `color_col`. If not
#' supplied, it is assumed to follow the order of the mean time of each roup in
#' `color_col`
#' @param abs_min starting point of the time line. If NULL, it is assumed to be
#' the minimum time in `x`
#' @importFrom magrittr %>%
break_time_line <- function (x, time_col, color_col, time_order=NULL,
                             abs_min=NULL, abs_max=NULL){
        x %>% dplyr::group_by (!!as.symbol (color_col) ) %>% 
              dplyr::summarise (min_pt=min(!!as.symbol(time_col)), 
                                mean_pt=mean(!!as.symbol(time_col)), 
                                max_pt=max(!!as.symbol(time_col))) %>%
              data.frame ()-> time_point
        if (is.null(time_order)) {
                 time_point %>% dplyr::arrange (mean_pt) %>%
                         dplyr::select (!!as.symbol (color_col)) %>%
                         tibble::deframe () -> time_order
        }

        N <- length(time_order)
        if (!is.null(abs_min)){ time_point$min_pt [time_point [, color_col] == time_order[1] ] <- abs_min}
        if (!is.null(abs_max)){ time_point$max_pt [time_point [, color_col] == time_order[N]] <- abs_max}
        # for each group, except for the first one whose min is certain for the
        # calculation above
        for (i in 2:N ){
                time_ind <- time_point [, color_col] %in% time_order[i]
                time_ind_p <- time_point [, color_col] %in% time_order[i-1]
                sel_time <- time_point$min_pt [time_ind]
                sel_time_p <- time_point$min_pt [time_ind_p]

                if (sel_time < sel_time_p){
                        print (paste (time_order[i], 'starts earlier than its previous group') )
                        time_point$min_pt [time_ind] <- time_point$max_pt [time_ind_p]
                }

                sel_time <- time_point$max_pt [time_ind]
                sel_time_p <- time_point$max_pt [time_ind_p]

                if (sel_time < sel_time_p){
                        print (paste (time_order[i-1], 'ends later than its next group') )
                        time_point$max_pt [time_ind_p] <- time_point$min_pt [time_ind]
                }
        }
        # resetting the maximum
        for (i in 2:N){
                time_ind <- time_point [, color_col] %in% time_order[i]
                time_ind_p <- time_point [, color_col] %in% time_order[i-1]
                time_point$max_pt [time_ind_p] <- time_point$min_pt [time_ind]
        }
        return (time_point %>% arrange (mean_pt))
}

#' Append the time interval labels
#'
#' @description This function uses the interval start and end points generated
#' from `break_time_line` to assign a time label for a vector
#' @param vec a vector to be labelled
#' @param time_df input dataframe to `break_time_line`
#' @param ... arguments to pass to `break_time_line`
append_break_labels <- function (vec, time_df, time_col, color_col,
                                 save_dir=NULL, save_label=NA,...){
        time_point <- break_time_line (time_df, time_col, color_col, 
                                       abs_min=min(vec), abs_max=max(vec))
        if (!is.null (save_dir) & !is.na (save_label)){
                file_name <- paste (save_dir, '/', time_col, '_', save_label, '.csv', sep='')
                utils::write.csv (time_point, file_name) 
        }
        lab_vec <- as.character (vec)
        for (i in 1:nrow(time_point)){
                lab_vec [vec > time_point$min_pt[i] & 
                         vec <= time_point$max_pt[i]] <- as.character (time_point [i, color_col])
        }
        lab_vec [vec == time_point$min_pt[1]] <- as.character (time_point [1, color_col])
        return (lab_vec)
}

#' Generate illustration plot for pseudotime progression
#'
#' @description This function is very specific to this paper and is not useful
#' for general applications.
#' @param cell_df a dataframe containing cell type information along the
#' pseudotime
#' @param time_df a dataframe containing names for pseudotime intervals
#' @param ... other arguments for `plot_two_branch`
#' @return a ggplot object
#' @export
graph_abs_time <- function (cell_df, time_df, AP=NULL, off_set=8, thickness=0.01,...){
        if (is.null(AP)){AP <- list(font_fam='sans', arrow_thickness=2)}
        plot_two_branch (cell_df, 
                         select_cells2 = cell_df$broad_type != 'STB',
                         select_cells1 = cell_df$broad_type != 'EVT',
                         text_offset=0.02, AP=AP, thickness=thickness,
                         off_set= off_set, time_axis=F, vert_trans=2*thickness,
                         time_df=time_df, time_time_col='peak_time', 
                         time_color_col='group', time_group_col='epil_branch',
                         branch_names=c('EVT_branch', 'STB_branch'), plot_branch_names=F,
                         repel_force=1,...
        ) -> time_plots
        ggplot2::ggplot () +
        time_plots [!1:length(time_plots) %in% c(2,4)]+
        plot_two_branch (cell_df, 
                         select_cells2 = cell_df$broad_type != 'STB',
                         select_cells1 = cell_df$broad_type != 'EVT',
                         text_offset=-0.005, AP=AP, thickness=thickness,
                         off_set= off_set, time_axis=T, 
                         branch_names=c('EVT_branch', 'STB_branch', 'TB_stem'),
                         ...
        )+ 
        time_plots [c(2,4)]+ #move the geom_text layers to the last
        # to prevent being overlapped by their plot elements
        ggplot2::theme (legend.position='none') 
}
