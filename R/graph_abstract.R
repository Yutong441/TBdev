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
sigmoid <- function (x, slope, off_set){
        return (1/( 1 +exp (-(x- off_set)/slope) ))
}

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
#' @importFrom ggplot2 aes_string
plot_one_branch <- function (dat, branch_id, branch_col='epil_branch', 
                             slope=1, off_set=5, 
                             filter_lim=c(-Inf, Inf),
                             time_col='MGP_PT', br_sign=1,
                             color_col='broad_type', thickness=0.01,
                             select_cells=NULL, label_group=NULL,
                             text_offset=0, AP=NULL){
        AP <- return_aes_param (AP)
        if (is.null(select_cells)){select_cells <- 1:nrow (dat)}
        dat <- dat [select_cells, ]
        dat [dat [, branch_col] %in% branch_id, ] %>%
                dplyr::filter (!!as.symbol (time_col) > filter_lim [1] ) %>%
                dplyr::filter (!!as.symbol (time_col) < filter_lim [2] ) %>%
                dplyr::mutate (y_ax = sigmoid (!!as.symbol (time_col), slope, off_set ) ) %>%
                dplyr::mutate (y_ax = y_ax * br_sign) %>% 
                dplyr::mutate (y_min = y_ax - thickness) %>%
                dplyr::mutate (y_max = y_ax + thickness) -> plot_data
        if (br_sign == 1){
                plot_data %>% dplyr::mutate (y_min = pmax (y_min, 0)) -> plot_data
        }else if (br_sign == -1){
                plot_data %>% dplyr::mutate (y_max = pmin (y_max, 0)) -> plot_data
        }

        if (!is.null(label_group)){ 
                dat %>% dplyr::filter (!!as.symbol (color_col) %in% label_group) -> lab_dat
        }else{lab_dat <- dat}
        lab_dat %>% dplyr::group_by (!!as.symbol (color_col) ) %>%
                dplyr::summarise (mean_pt = mean (!!as.symbol (time_col) )) %>%
                dplyr::mutate (y_max = sigmoid (mean_pt, slope, off_set )+ 
                        thickness + text_offset) %>%
                dplyr::mutate (y_max = y_max*br_sign) -> label_dat
        layer1 <- ggplot2::geom_ribbon (aes_string (x=time_col, ymin='y_min', ymax='y_max', 
                                            fill=color_col), data=plot_data) 
        layer2 <- ggplot2::geom_text ( aes_string (x='mean_pt', y='y_max', label=color_col), 
                             data = label_dat, inherit.aes=F,
                             size=AP$point_fontsize, family=AP$font_fam)
        return (list (layer1, layer2))
}

#' Draw 2 sigmoid curves
#' 
#' @param branch1 which samples in `dat` make up the first branch
#' @param select_cells1 which cells in the first branch
#' @param time_axis whether to plot the x axis
#' @param axis_label what to label the `time_axis`
#' @seealso `plot_one_branch`
#' @examples
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
                             branch1 = c('main', 'EVT_b'), 
                             branch2= c('main', 'STB_b'), 
                             select_cells1 = NULL,
                             select_cells2 = NULL, AP=NULL,
                             time_axis=T, axis_label='time',...){
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

        ggplot2::ggplot () +
                plot_one_branch (dat, branch2, time_col = time_col, 
                                 color_col = color_col, select_cells= select_cells2,
                                 branch_col=branch_col, br_sign=-1, 
                                 filter_lim=filter_lim, label_group=label_group, AP=AP,...)+
                plot_one_branch (dat, branch1, time_col = time_col, 
                                 color_col = color_col, select_cells= select_cells1,
                                 branch_col=branch_col, filter_lim=filter_lim, AP=AP,...)+
                theme_TB ('no_arrow', feature_vec=dat[, color_col], color_fill=T, AP=AP) -> plot_ob
        if (time_axis){
                plot_ob +ggplot2::geom_segment (aes (x=min_b, y=0, 
                         xend=max_b*0.85, yend=0), arrow=get_arrow (AP),
                         size=AP$arrow_thickness, linejoin = AP$arrow_linejoin)+
                         ggplot2::annotate('text', x=max_b*0.9, y=0, label=axis_label, 
                                       size=AP$point_fontsize, family=AP$font_fam, 
                                       hjust='left') -> plot_ob
        }
        return (plot_ob)
}
