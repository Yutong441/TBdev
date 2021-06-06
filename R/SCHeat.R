# heatmap for scRNA-seq

get_rainbow_col <- function (named_vector, AP, default=T, provided_color=NULL,
                             regexp=NULL){
        all_matches <- unique (named_vector) %in% names (provided_color)
        if (!is.null (provided_color) & mean (all_matches) == 1 ){
                print ('entering provided color')
                col_color <- provided_color
        }else{
                if (default){
                        label_group <- gtools::mixedsort (unique (named_vector) )
                        col_color <- grDevices::rainbow (length (label_group))
                        names (col_color) <- as.character(label_group)
                        names (col_color) [is.na (names (col_color))] <- 'NA'
                }else{
                        col_color <- custom_color (named_vector, AP,regexp=regexp)
                        print (col_color)
                }
        }
        return (col_color)
}

#' @importFrom grid gpar unit
#' @noRd
get_heat_param <- function (AP, title_pos, grid_height){
        list (labels_gp=gpar(fontsize=AP$fontsize, fontfamily=AP$gfont_fam),
              title_gp=gpar(fontsize=AP$fontsize, fontfamily=AP$gfont_fam, fontface='italic'),
              title_position=title_pos, 
              grid_height= unit (grid_height, 'mm')
        )
}

get_hori_bars <- function (x, group.by, group_order, default_color,
                           provided_color, reorder_column,
                           column_reorder_levels, AP){
        color_map_list <- list ()
        hori_bar <- list ()
        for (i in 1:length (group.by)){
                hori_data <- x@meta.data[, group.by[i] ]
                hori_bar[[i]] <- as.vector (hori_data[group_order])
                column_color <- get_rainbow_col (hori_data[group_order],AP,
                                                 default=default_color,
                                                 provided_color=provided_color)
                color_map_list [[i]] <- column_color
        }

        print ('ordering columns')
        names (color_map_list) <- group.by
        HA_df <- data.frame (do.call(cbind, hori_bar ))
        colnames (HA_df) <- group.by
        rownames (HA_df) <- colnames (x) [group_order]
        if (reorder_column){
                for (i in group.by){
                        HA_df [,i] <- partial_relevel (HA_df [,i], column_reorder_levels)
        }}
        return (list (hori_bar, HA_df, color_map_list) )
}

return_row_HA_ob <- function (df_list, col_list, vert_anna_param, AP, show_legend=T){
        ComplexHeatmap::HeatmapAnnotation (df=df_list, which='row', col=col_list,
                           show_annotation_name=F,
                           annotation_legend_param=vert_anna_param,
                           show_legend=show_legend,
                           annotation_name_gp = grid::gpar (fontsize=AP$fontsize, 
                                                      fontfamily=AP$gfont_fam)
        ) 
}

#' Perform row or column scaling
#'
#' @param xx a numerical matrix or dataframe
#' @importFrom magrittr %>%
scaling <- function (xx, row_scale, column_scale){
        xx <- as.matrix ( xx)
        if (row_scale){ xx %>% t() %>% scale () %>% t() -> xx }
        if (column_scale){ xx %>% scale () -> xx }
        xx [is.na (xx) ]<- 0
        return (xx)
}

#' @export
scale_seurat <- function (x, row_scale=F, column_scale=F, slot_data='data',
                          assay='RNA'){
        exp_mat <- Seurat::GetAssayData (x, slot=slot_data, assay=assay)
        exp_mat <- scaling (exp_mat, row_scale, column_scale)
        return (Seurat::SetAssayData (x, slot=slot_data, assay=assay, new.data=exp_mat))
}

#' @export
scale_seurat_with_seurat <- function (x, ref, slot_data=c('data', 'data'), 
                                      assay=c('RNA', 'RNA'), row_scale=F, 
                                      column_scale=F){
        exp_mat <- as.matrix (Seurat::GetAssayData (x, slot=slot_data[1], assay=assay[1]))
        exp_ref <- as.matrix (Seurat::GetAssayData (ref, slot=slot_data[2], assay=assay[2]))
        if (row_scale){
                mean_ref <- rowMeans (exp_ref)
                sd_ref <- apply (exp_ref, 1, stats::sd)
                exp_mat <- (exp_mat - mean_ref )/sd_ref
                exp_mat [is.na (exp_mat)] <- NA
        }
        if (column_scale){
                print ('column scaling is not supported')
        }
        return (Seurat::SetAssayData (x, slot=slot_data[1], assay=assay[1], new.data=exp_mat))
}

#' Obtain color gradient for heatmap
#'
#' @param xx a matrix or dataframe that will be the input of the heatmap
#' @param quantile_val the dynamic range of the heatmap. For example, if it is
#' 0.05, the heatmap gradient will go from 5% to 95% percentile of the values
#' in `xx`
#' @param center_scale whether to color the middle of the dynamic range as
#' white. Otherwise, 0 is colored white.
#'
#' @importFrom stats quantile
#' @importFrom circlize colorRamp2
determine_color_gradient <- function (xx, quantile_val=0., center_scale=T,
                                      AP=NULL){
        AP <- return_aes_param (AP)
        color_vec <- AP$heatmap_color
        plot_min <- quantile(xx, quantile_val)
        plot_max <- quantile(xx, 1-quantile_val)

        if (center_scale){center <- (plot_min+plot_max)/2}else{center <- 0}
        if (min (xx) < 0){
                if (center_scale){
                        plot_min <- floor (plot_min)
                        plot_max <- ceiling(plot_max)
                        peak_mag <- abs (max(-plot_min, plot_max))
                        plot_min <- -peak_mag
                        plot_max <- peak_mag
                        center <- 0
                }
                color_scale<- colorRamp2 (c(plot_min, center, plot_max ), color_vec)
        }else{
                color_scale<- colorRamp2 (round (c(plot_min, center, 
                                            plot_max), 1), color_vec)
        }
        break_points <- round (c(plot_min, center, plot_max ),1)
        return (list (color_scale, break_points))
}

make_heat_legend <- function (anna_param, col_fun, name=NA){
        main_anna_param <- append (anna_param, list (title=name,
                                                     col_fun=col_fun) )
        return (do.call (ComplexHeatmap::Legend, main_anna_param))
}

#' Create annotation legend manually
#'
#' @description This function takes in similar arguments as `HeatmapAnnotation`
#' and output the legend only. This is to give a better graphics control over
#' the legend position.
#' @param anna_param annotation param for `annotation_legend_param`
#' @param legend_bar_title what to appear in the title of legend bars in the
#' order of appearance
#' @param color the color for those titles in `legend_bar_title`
#' @param name name of the legend
make_anno_legend <- function (anna_param, legend_bar_title, color, name=NA){
        row_fill_color <- color [match (legend_bar_title, names (color) ) ]
        row_param <- append (anna_param, list (at=legend_bar_title, 
                                               legend_gp=gpar(fill=row_fill_color),
                                               title=name) )
        return (do.call (ComplexHeatmap::Legend, row_param))
}

make_multi_anno_legend <- function (anna_param, anno_df, color_map_list, anno_names, index){
        make_anno_legend (anna_param, levels (anno_df [, index] ),  
                          color_map_list[[index]], anno_names[index])
}

#' Heatmap for Seurat object
#'
#' @description Equivalent to `Seurat::DoHeatmap` with the addition of sidebar
#' along the rows (features) to group them using ComplexHeatmap. However,
#' currently, only one sidebar can be plotted for rows. There is no limit for
#' the number of side bars for columns.
#' 
#' @param x Seurat object
#' @param group.by by which feature are the cells grouped
#' @param color_row a vector of gene names with the name of the vector elements
#' being the lineages
#' @param assay the Seurat assay in `x`
#' @param slow_data the slot in `assay`
#'
#' @param highlight which row names to highlght
#' @param highlight_names display names for the highlighted and non-highlighted
#' groups
#' @param show_column_names whether to display column names (Don't confuse it
#' with column title, which is group name)
#' @param show_row_names whether to show row names
#' @param show_column_anna whether to show the category for the group names for
#' the columns around the column color bar. There is not a similar option for
#' rows because this function only supports one row sidebar
#' @param column_rotation rotate column name labels, 90 is vertical, 0 is
#' horizontal, there are no other choices
#' @param annotation_name_side side of the annotation name for column side bar
#' @param show_column_bars whether to show horizontal bars for the columns. It
#' can be a vector of the same length as `group.by` to specify the group(s) to
#' show the horizontal bars.
#'
#' @param cluster_columns whether to perform hierarchical clustering on columns
#' @param cluster_rows whether to perform hierarchical clustering on rows
#' @param column_split whether to split the columns belonging to the same group
#' in `group.by` by a thin dividing line. '1' indicate splitting. 'NA'
#' indicates no spliting
#' 
#' @param row_reorder_levels a vector indicating the order of labels for row
#' group names that would appear from top to bottom of the heatmap. The default
#' is `cell_order` from the list supplied to the `AP` argument
#' @param column_reorder_levels a list of vectors. Each item in the list
#' indicates the order of each column color bar labels. Again the default is in
#' the `AP` argument
#' @param reorder_columns whether to order the columns according to
#' `column_reorder_levels`
#'
#' @param row_scale whether to perform row scaling.
#' @param column_scale whether to perform column scaling.
#' @param center_scale whether to put white colors half way between the min and
#' max of the matrix. Otherwise, white color is placed at value 0.
#' @param color_scale the scale of color gradient for heatmap. By default, the
#' `heatmap_color` field of the list supplied to `AP` should define the color
#' for min, middle and max of the dataset. Otherwise, please input a function
#' created by `circlize::colorRamp2`
#' @param default_color whether to use the default color settings of
#' ComplexHeatmap for main matrix
#' @param quantile_val if supplied, then the color gradient would not extend
#' over the entire range of the matrix. The default 0.05 means that color
#' gradient would extend from 5% to 95% percentile. This is to avoid outliers
#' skewing the dynamic range of colors.
#' @param provided_color a named vector where the vector contains the colors
#' and the names indicate the group labels that would be colored by the
#' corresponding color in the sidebar.
#' @param AP a list of aesthetic parameters. The important fields include
#' `cell_order` for ordering group labels, `gfont_fam`, `fontsize` and
#' `heatmap_color`
#' 
#' @param title_pos adjustment of legend title
#' @param heat_name name of the main heatmap to appear in the legend
#' @param break_points where to show values in the heatmap gradient legend
#' @param row_legend_labels name of the row sidebar to appear in the legend
#' @param column_legend_labels name of the column sidebar in the legend. You may
#' supply as many items as the number of column sidebars
#' @param left_HA whether to show row sidebars
#' @param top_HA whether to show column sidebars
#' @param row_titles whether to show row titles. choose NULL if not
#' @param group_order the ordering of each column
#' @param main_width width of the heatmap in cm
#' @param main_height height of the heatmap in cm
#' @param grid_height legend grid height, increase this value if the spacing
#' between legend labels become too narrow
#' @param heat_grid_height grid height for the color gradient legend. If NULL,
#' it is the same as `grid_height`
#' @param automatic whether the default legend generation process in
#' ComplexHeatmap is to be used. If not, all legends will be aligned vertically
#' and only be arranged in separate columns if the total length exceeds the
#' heatmap length. 
#'
#' @importFrom ComplexHeatmap HeatmapAnnotation Heatmap Legend
#' @importFrom grid gpar unit
#' @importFrom stats quantile
#' @author Yutong Chen
#' @export
seurat_heat <- function (x, group.by=NULL, color_row=NULL, 
                         # Seurat related
                         assay='RNA', slot_data='data', 

                         # displaying row and column names
                         highlight=NULL, highlight_names=NULL, 
                         show_column_names=F, show_row_names=T,
                         show_column_anna=T, column_rotation=0,
                         annotation_name_side='left',
                         column_names_side = 'bottom',
                         column_title_side = 'top',
                         row_names_side = 'left',
                         row_title_side = 'left',
                         show_column_bars = T,
                         row_title_fontface='bold',
                         column_title_fontface='bold',

                         # clustering
                         cluster_columns=F,
                         cluster_rows=F, 
                         column_split=1, 

                         # ordering row and column name
                         row_reorder_levels=NULL, column_reorder_levels=list(),
                         reorder_column=T, 

                         # color setting
                         row_scale=F, column_scale=F,
                         center_scale=F, 
                         color_scale=NULL,
                         default_color=F, 
                         quantile_val=0.05,
                         provided_color=NULL, 
                         AP=NULL, 

                         # legend
                         title_pos='topleft',
                         heat_name='norm count', break_points=NULL,
                         row_legend_labels='DE genes',
                         column_legend_labels=NULL, 
                         column_legend_order=NULL,
                         left_HA=T, top_HA=T,
                         row_titles=character(0),
                         group_order=NULL, 
                         main_width=NULL, main_height=NULL,
                         grid_height=8,
                         heat_grid_height=NULL,
                         automatic=T,
                         ...){

        AP <- return_aes_param (AP)
        if (is.null(color_row)){color_row <- rownames (x)}
        # disable all legends if automatic option is FALSE
        if (!automatic){
                show_column_legend <- F
                show_row_legend <- F
                show_heat_legend <- F
        }else{
                show_column_legend <- T
                show_row_legend <- T
                show_heat_legend <- T
        }

        if (is.null (group.by)){
                x$any_group <- 'none'
                group.by <- 'any_group'
                column_split <- NA
                top_HA <- F
        }

        # setting default row and column order levels
        if (is.null (column_legend_labels)){column_legend_labels <- group.by}
        if (is.null(row_reorder_levels) ){row_reorder_levels <- AP$cell_order}
        if (length(column_reorder_levels) == 0){column_reorder_levels <- list (AP$cell_order)}
        inp_col <- length (column_reorder_levels)
        inp_group <- length (group.by)

        if (inp_col != inp_group){
                print ('adding column sorting levels')
                added_lev<- rep (column_reorder_levels [inp_col], inp_group - inp_col) 
                column_reorder_levels <- c(column_reorder_levels, added_lev)
        }
        column_reorder_levels <- unlist (unique (column_reorder_levels) )

        # check rownames
        unique_index <- !duplicated (color_row) & color_row %in% rownames (x) 
        color_row <- color_row [unique_index]
        if (is.null (names (color_row )) ){
                names (color_row) <- rep ('NA', length (color_row)) 
        }
        color_row_names <- names (color_row)
        color_row_names [is.na (color_row_names)] <- 'NA'
        color_row_names <- partial_relevel (color_row_names, row_reorder_levels)
        if (mean (color_row_names == 'NA') == 1){
                color_row_names <- rep ('no values', length(color_row_names))
                left_HA <- F
                row_titles <-NULL 
        }

        # order vectors
        if (reorder_column){
                for (i in 1:length(group.by)){
                        one_feature <- x@meta.data[, group.by[i] ]
                        x@meta.data[, group.by[i]] <- partial_relevel (x@meta.data[, group.by[i]], 
                                                                     column_reorder_levels)
        }}
        if (length (group.by) > 1){
                x@meta.data %>% dplyr::select (group.by) %>%
                        tidyr::unite ('joined', group.by) -> indiv_names
                indiv_names <- indiv_names$joined
        }else{indiv_names <- x@meta.data [, group.by]}
        if (reorder_column) {indiv_names <- partial_relevel (indiv_names, column_reorder_levels)}
        if (is.null(group_order)) {group_order <- order (indiv_names)}

        print ('making horizontal bar')
        HA_all_list <- get_hori_bars (x, group.by, group_order, default_color,
                                   provided_color, reorder_column,
                                   column_reorder_levels, AP)
        hori_bar <- HA_all_list[[1]]
        HA_df <- HA_all_list[[2]]
        color_map_list <- HA_all_list[[3]]

        print ('legend parameter settings')
        anna_param <- get_heat_param (AP, title_pos, grid_height)
        if (!is.null(column_legend_labels)){
                colnames (HA_df) <- column_legend_labels
                names (color_map_list) <- column_legend_labels
        }
        print ('finally making horizontal bars')
        # NB: HeatmapAnnotation may fail system check when in `do.call`
        if (length(show_column_bars)==1){
                show_column_bars <- rep (show_column_bars, length(group.by))
        }
        hori_bars<- ComplexHeatmap::columnAnnotation (df=HA_df [, show_column_bars, drop=F], 
                                       col= color_map_list [show_column_bars],
                                       show_annotation_name=show_column_anna,
                                       show_legend =show_column_legend,
                                       annotation_name_side=annotation_name_side,
                                       annotation_legend_param=anna_param,
                                       annotation_name_gp = gpar (fontsize=AP$fontsize, 
                                                                  fontfamily=AP$gfont_fam)
        )

        print ('making vertical bar')
        vert_anna_param <- append (anna_param, list (title=row_legend_labels))
        if ( is.null (highlight) | is.null (highlight_names)  ){
                row_color <- get_rainbow_col (color_row_names ,AP,
                                              default=default_color, regexp='\n(.*)$')

                vert_bars <- return_row_HA_ob (data.frame (Marker = color_row_names), 
                                             list(Marker=row_color),
                                             vert_anna_param, AP,
                                             show_legend=show_row_legend)
        }else{
                row_color <- get_rainbow_col (color_row_names, AP,
                                              default=default_color)
                highlight <- rev (highlight_names) [as.factor (highlight)] [unique_index]
                label_high <- c('blue', 'green')
                names (label_high) <- highlight_names
                vert_bars <- return_row_HA_ob (data.frame (Marker = color_row_names,
                                                         highlight = highlight),
                                              list(Marker=row_color, highlight=label_high),
                                              vert_anna_param, AP, show_legend=F) 
        }

        print ('processing expression matrix')
        plot_data <- Seurat::GetAssayData (x, slot=slot_data, assay=assay)
        plot_data <- plot_data [match (color_row, rownames (plot_data)), group_order]
        plot_data <- scaling (plot_data, row_scale, column_scale)

        print ('setting heatmap controls')
        if (!is.na(column_split)){column_split <- HA_df [, column_split] 
        }else{column_split <- NULL}
        if (!is.null(main_width)){main_width <- unit (main_width, 'cm') }
        if (!is.null(main_height)){main_height <- unit (main_height, 'cm') }
        if (is.null(color_scale)){
                color_grad <- determine_color_gradient (plot_data, quantile_val, 
                                                        center_scale, AP)
                color_scale <- color_grad [[1]]
                break_points <- color_grad [[2]]
                main_legend_param <- append(anna_param, list (at=break_points))
        }else{
                if (!is.null (break_points) ) {
                        main_legend_param <- append (anna_param, list (at=break_points) )
                }else{main_legend_param <- anna_param}
        }
        if (is.null(heat_grid_height)){heat_grid_height <- grid_height}
        main_legend_param$grid_height <- unit (heat_grid_height, 'mm')

        if (!top_HA){hori_bars <- NULL}
        if (!left_HA){vert_bars <- NULL}

        print ('drawing heatmap')
        Heatmap (plot_data, cluster_rows=cluster_rows, 
                 cluster_columns=cluster_columns, 
                 show_column_names=show_column_names,
                 column_names_side=column_names_side,
                 column_title_side=column_title_side,
                 column_names_gp=gpar (fontsize=AP$fontsize, 
                                       fontfamily=AP$gfont_fam),
                 show_row_names=show_row_names,
                 top_annotation=hori_bars, 
                 left_annotation=vert_bars,

                 column_split=column_split,
                 row_split=color_row_names, 
                 row_names_side=row_names_side, 
                 row_title_side=row_title_side,
                 row_names_gp = gpar (fontsize=AP$fontsize, 
                                      fontfamily=AP$gfont_fam), 

                 row_title=row_titles,
                 row_title_rot=0,
                 column_title_rot=column_rotation,
                 row_title_gp = gpar (fontsize=AP$fontsize, fontfamily=AP$gfont_fam, 
                                      fontface=row_title_fontface),
                 column_title_gp = gpar (fontsize=AP$fontsize, fontfamily=AP$gfont_fam, 
                                         fontface=column_title_fontface), 

                 col=color_scale,
                 width = main_width, height = main_height,
                 name=heat_name,
                 heatmap_legend_param = main_legend_param, 
                 show_heatmap_legend = show_heat_legend, ...
                 )  -> heat_ob

        if (!automatic){
                print ('drawing legends')
                heat_legend <- make_heat_legend (main_legend_param, color_scale, heat_name)

                if (left_HA){
                        row_leg_lab <- gsub ('\n(.*)$', '', levels(color_row_names))
                        names (row_color) <- gsub ('\n(.*)$', '', names(row_color))
                        row_legend <- make_anno_legend (anna_param, row_leg_lab, 
                                                             row_color, row_legend_labels )
                }else{row_legend <- NULL}
                if (top_HA ){
                        column_legends <- lapply (as.list (1:length (column_legend_labels)), 
                                                  function (x){make_multi_anno_legend (anna_param, 
                                                HA_df, color_map_list, column_legend_labels, x) } )
                }else{column_legends <- list (NULL) }
                if (is.null(column_legend_order)){
                        column_legend_order <- 1:length (show_column_bars)
                }
                legend_list<- append (list(row_legend, heat_legend), column_legends [show_column_bars] [column_legend_order])
                legend_list<- legend_list [sapply (legend_list, function(x){!is.null(x)})]
                pd <- ComplexHeatmap::packLegend (list=legend_list,
                                                  max_height=main_height,
                                                  direction='vertical')
                ComplexHeatmap::draw (heat_ob, annotation_legend_list= pd)
        }else{
                return (heat_ob)
        }
}
