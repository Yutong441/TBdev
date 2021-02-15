# This script contains the functions to make a poster
# If the poster contains many vector graphics and very few texts, this script is ideal
# because loading a lot of vector graphics into inkscape or adobe illustrator
# can slow their program extensively.
# The prinicple behind coding a poster in R is by breaking down the page into a
# grid of cells. Then the content for each cell is specificied.
#-------------------------------------------------- 

# ----------Managing text style----------
get_default <- function (){
        list (fontface='plain', fontsize=35, font_fam='arial', xpos=0,
              ypos=0.95, hjust=0, style='default')
}

stylize <- function (style_name){
        if (style_name=='default'){get_default()
        }else if (style_name == 'heading') {
                list (fontface='bold', fontsize=48, font_fam='arial', xpos=0,
                      ypos=0.95, hjust=0, style='heading')
        }else if (style_name == 'title'){
                list (fontface='bold', fontsize=70, font_fam='arial', xpos=0.5,
                      ypos=0.5, hjust='center', style='title')
        }else if (style_name == 'subheading'){
                list (fontface='italic', fontsize=38, font_fam='arial', xpos=0,
                      ypos=0.95, hjust=0, style='subheading')
        }
}

#' Customise text appearance from `grid.text`
#'
#' @param text_str a single string
#' @param xpos position in the text horizontally as a percentage of the width
#' @param ... arguments for `grid::grid.text`
custom_text <- function (text_str, fontface, fontsize, font_fam, xpos, ypos, ...){
        cus_font <- get_highlight_font (fontface, fontsize, font_fam)
        grid::grid.text (text_str, 
                         x=grid::unit (xpos, 'npc'),
                         y=grid::unit (ypos, 'npc'), 
                         gp=cus_font, ...)
}

#' Customise text using an argument string
#'
#' @param text_args a string of keyword argument separated by comma, e.g.
#' 'fontsize=5, fontfam='arial''
#' @importFrom magrittr %>%
custom_text_str <- function (text_str, text_args){
        # convert the `text_args` into list
        text_args <- paste ('list(', text_args, ')')
        parse(text=text_args) %>% eval() -> arg_list
        arg_list [['text_str']] <- text_str

        # append default parameters
        if (!'style' %in% names (arg_list) ){arg_list$style <- 'default'}
        style_list <- stylize (arg_list$style)
        arg_list <- return_aes_param (arg_list, style_list)
        arg_list <- arg_list [names (arg_list) !='style']
        do.call (custom_text, arg_list)
}

# ----------configuration array----------

#' Obtain row start and end coordinates
row_start_end_coord <- function (config_arr){
        config_arr$row.start <- 1
        config_arr %>% dplyr::filter (block==1 & row.order==1) %>% 
                dplyr::slice_max (row.num, n=1) %>% dplyr::pull (row.num) -> abs_len
        for (i in unique (config_arr$block)){
                # in a particular block
                sel_index <- config_arr$block==i
                sel_arr <- config_arr [sel_index, ]
                for (j in 1:nrow (sel_arr)){
                        ind <- sel_arr$row.order[j] 
                        # Do not bother with the first row
                        if (ind != 1){
                                prev_ind <- which (sel_arr$row.order==ind-1)
                                if (length(prev_ind)!=0){
                                        prev_start <- sel_arr$row.start [prev_ind]
                                        prev_len <- sel_arr$row.num [prev_ind]
                                        sel_arr$row.start [j] <- max(prev_start + prev_len)
                                }else{
                                        print ('cannot find the previous plot')
                                        sel_arr$row.start [j] <- 1+ abs_len
                                }
                        }
                }
                config_arr$row.start [sel_index] <- sel_arr$row.start 
        }
        config_arr$row.end <- config_arr$row.start + config_arr$row.num-1
        return (config_arr)
}

base_layout <- function (config_arr){
        M <- max (config_arr$row.end) - min (config_arr$row.start)+1
        N <- max (config_arr$column.end) - min (config_arr$column.start)+1
        return (matrix (0, M, N))
}

create_layout <- function (grob_list, config_arr){
        base_arr <- base_layout (config_arr)
        no_coord <- names (grob_list) [!names (grob_list) %in% config_arr$label]
        if (length (no_coord)!=0){
                warning (paste ('The following objects lack coordinates \n', no_coord))
        }
        for (i in 1:length (grob_list)){
                ID <- names (grob_list)[i]
                if (!ID %in% config_arr$label){
                        stop (paste ('No coordinates found for', ID, 
                                     '. Check config_arr'))
                }
                x1 <- config_arr$row.start [config_arr$label==ID]
                x2 <- config_arr$row.end [config_arr$label==ID]
                y1 <- config_arr$column.start [config_arr$label==ID]
                y2 <- config_arr$column.end [config_arr$label==ID]
                base_arr [x1:x2, y1:y2] <- i
        }
        return (base_arr)
}

# ----------putting everything together----------

#' Put text and figures into a poster
#'
#' @param grob_list a named list of objects to plot
#' @param save_path where to save the poster
#' @param config_arr a data.frame with the following columns:
#' `label`: label for the item that matches the name of the item in `grob_list`
#' `row.order`: order of the item along the row
#' `row.num`: number of rows
#' `column.start`: from which column the item span starts
#' `column.num`: number of columns
#' `text.style`: a string specifying the style of the text, available options
#' are 'fontsize', 'font_fam', 'fontface', 'xpos', 'ypos', 'hjust', 'vjust' and
#' any arguments that pass into `grid::grid.text`
#'
#' The function calculate internally: 
#' `column.end`: from which column the item span ends
#' `row.start`: from which row the item span starts
#' `row.end`: from which row the item span ends
#' @description Empirally, I find one figure best correspond to 12 lines
#' @return a saved pdf file
#' @importFrom grid pushViewport viewport grid.layout
#' @export
arrange_poster <- function (grob_list, save_path, config_arr, plot_width=7,
                            plot_height=7, margin_width=0.79,
                            margin_height=0.79, page_width=NULL,
                            page_height=NULL){
        config_arr <- row_start_end_coord (config_arr)
        config_arr$column.end <- config_arr$column.start+config_arr$column.num-1
        grid_layout <- create_layout (grob_list, config_arr)
        plot_width <- expand_length (plot_width, ncol (grid_layout)) # a vector
        plot_height <- expand_length (plot_height, nrow (grid_layout))
        if (is.null(page_width)){page_width <-  sum (plot_width) }
        if (is.null(page_height)){page_height <- sum (plot_height)}
        width_pro <- 1 - margin_width/page_width
        height_pro <- 1 - margin_height/page_height

        grDevices::cairo_pdf (save_path, width=page_width, height=page_height)
        grid::grid.newpage ()
        pushViewport(viewport(layout=grid.layout(nrow (grid_layout), ncol(grid_layout)) ,
                              width=width_pro, height=height_pro))

        # make sure every object in `grob_list` has been assigned a
        # corresponding spot in `grid_layout`
        Len <- length (grob_list)
        grid_ind <- unique (c(grid_layout))
        no_grid_name <- names (grob_list)[!1:Len %in% grid_ind]
        if (length (no_grid_name) == 0){
                warning (paste ('The following objects have no specified grids \n', no_grid_name))
        }
        # start plotting
        for (i in 1:Len){
                row_pos <- unique (which (grid_layout == i, arr.ind=T) [, 'row'])
                col_pos <- unique (which (grid_layout == i, arr.ind=T) [, 'col'])
                print (paste ('arranging figure', i, names (grob_list)[i] ))
                pushViewport (viewport (layout.pos.col=col_pos,layout.pos.row=row_pos))
                text_arg <- config_arr$text.style[config_arr$label == names (grob_list)[i]]
                push_plot_viewport (grob_list[[i]], text_arg)
                popViewport (1)
        }
        grDevices::dev.off ()
}
