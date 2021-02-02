#' find the location of the maximum value
#'
#' @param x a dataframe with time series information along the columns, and the
#' features along the rows. It must contain one row that has the time index.
#' Otherwise, the time index is initiated as from 1 to the sample number
#' @param max_only detect absolute peak, instead of distinguishing between
#' positive or negative peak
#' @return a dataframe with the following columns:
#' `val`: magnitude of the peak relative to the trough
#' `feature`: genes
#' `change_sign`: whether it is a positive peak or negative peak
#' `peak_time`: when the peak occurs
#' @importFrom magrittr %>%
#' @export
find_peak <- function (x, time_index=NULL, max_only=T){
        x %>% t () %>% as.data.frame () -> x_df
        if (is.null (time_index)){
                time_index <- 'time_ind'
                x_df  %>% tibble::add_column (time_ind =1:ncol (x_mat)) -> x_df
        }
        x_df %>% dplyr::select (!all_of (time_index) ) -> no_time 
        if (!max_only){no_time_abs <- abs (no_time) 
        }else{no_time_abs <- no_time}

        no_time_abs %>% t()  %>% max.col () -> peak_loc
        peak_loc <- x_df [peak_loc, time_index]
        peak_df <- data.frame ('feature'=colnames (no_time), 'peak_time'=peak_loc)
        # get Z score
        max_val <- no_time_abs %>% apply (2, max)
        min_val <- no_time_abs %>% apply (2, min)
        peak_df$val <- max_val - min_val
        rownames (peak_df) <- peak_df$feature

        if (!max_only){
                true_max <- no_time %>% apply (2, max)
                peak_df$change_sign <- max_val == true_max
        }
        return (peak_df)
}

#' Detect genes that show the peak changes in GP
#'
#' @param seurat_ob a Seurat object. One row must contain the information of
#' pseudotime for each cell, as in the `time_ind` argument.
#'
#' @export
find_transition_from_seurat <- function (seurat_ob, assay='RNA', slot_data='counts',
                                   time_ind='pseudotime'){
        exp_mat_inf1 <- as.matrix(Seurat::GetAssayData (seurat_ob, assay=assay, slot=slot_data) )
        exp_mat_diff <- apply (exp_mat_inf1, 1, diff)
        dpt <- (max(exp_mat_inf1['pseudotime',]) - min (exp_mat_inf1['pseudotime',]))/ncol(exp_mat_inf1)
        exp_mat_diff <- exp_mat_diff/dpt
        exp_mat_diff [, 'pseudotime'] <- exp_mat_inf1 ['pseudotime', 1:ncol(exp_mat_inf1)-1]
        return ( find_peak (t(exp_mat_diff), time_index='pseudotime', max_only=F))
}

#' @importFrom magrittr %>%
#' @importFrom ggplot2 aes_string
#' @noRd
plot_peak_line <- function (peak_df, AP, gene_col='feature', time_col='peak_time') {
        plot_dat <- peak_df
        colnames (plot_dat) [colnames (plot_dat) == gene_col ] <- 'gene'
        plot_dat$time_lab <- format (round (plot_dat [, time_col], 2), nsmall=2)
        layer1 <- ggplot2::geom_vline (aes_string (xintercept=time_col), data =
                                       plot_dat, linetype='dashed', size=1) 
        layer2 <- ggplot2::geom_text (aes_string (x=time_col, y=-Inf, label='time_lab'), 
                                      vjust=0.5, size=AP$point_fontsize, data=plot_dat)
        layer3 <- ggplot2::coord_cartesian (clip='off')
        return (list (layer1, layer2, layer3))
}

#' Save peak genes in a particular KEGG pathway
#'
#' @importFrom magrittr %>%
#' @export
save_peak_genes <- function (xx, save_path, pathway=NULL, pathway_df=NULL,
                             max_num=NULL, arrange_val='val'){
        if (!is.null (pathway) & !is.null (pathway_df)){
                sel_genes <- gene_per_term (pathway_df, pathway, return_val=T)
                xx %>% dplyr::filter (feature %in% sel_genes[[1]]) -> xx
        }
        if (!is.null (max_num)){
                xx %>% dplyr::slice_max (abs (!!as.symbol (arrange_val) ), n=max_num) -> x
        }
        xx %>% dplyr::arrange (dplyr::desc (abs (!!as.symbol (arrange_val) ) ) ) %>% 
                utils::write.csv (save_path)
}

#' @importFrom magrittr %>%
get_gene_one_group <- function (xx_df, group_name, gene_num){
        xx_df %>% dplyr::filter (group == group_name) %>% 
                dplyr::slice_max (val, n=gene_num)
}

#' @importFrom magrittr %>%
get_gene_all_groups <- function (xx_df, peak_df, min_genes){
        peak_df %>% dplyr::count (group) %>% tibble::deframe () -> count_group
        count_group <- round (count_group*min_genes/min(count_group))
        group_uniq <- as.list (unique (xx_df$group))
        xx_df_sel <- lapply (group_uniq, function (x){
                                get_gene_one_group (xx_df, x, count_group[x]) 
                             })
        return (do.call (rbind, xx_df_sel))
}

#' Get genes for labelling in heatmap
#'
#' @param xx a dataframe generated from `find_peak`. For the columns that it
#' should contain, see `find_peak`
#' @param plot_quantile the top quantile of genes in `xx` to show in heatmap
#' @param min_genes minimum number of genes to show in the heatmap row labels
#' @importFrom magrittr %>%
#' @export
peak_gene_for_heatmap <- function (xx, plot_quantile=0.95, min_genes=2){
        # get the genes to plot in heatmap
        xx %>% dplyr::filter (val > stats::quantile (xx$val, plot_quantile) ) %>% 
                arrange (desc(val)) -> label_peak

        label_peak %>% dplyr::arrange (change_sign*val ) %>% 
                dplyr::select (group, feature) %>% 
                dplyr::arrange (group) %>% 
                tibble::deframe () -> sel_genes

        # get the genes to show in heatmap row labels
        xx_df <- get_gene_all_groups (xx, label_peak, min_genes)
        xx_df %>% dplyr::group_by (group) %>%  
                dplyr::select (group, feature) %>%
                dplyr::summarise (all_genes = paste (feature, collapse=', ')) %>% 
                data.frame ()-> gene_list

        gene_list$all_genes <- paste (gene_list$group, '\n(', gene_list$all_genes, ')', sep='')
        gene_list$all_genes <- sapply (gene_list$all_genes, function(x){
                                        paste (line_break (strsplit (x, ',')[[1]]), collapse=',')
                             }) # see `line_break` in `enrichment.R`
        # remove the trailing comma
        gene_list$all_genes <- gsub (',$', '', gene_list$all_genes)
        gene_list$all_genes <- gsub (',\n$', '', gene_list$all_genes)
        names (sel_genes) <- gene_list$all_genes [match (names(sel_genes), gene_list$group)]
        return (sel_genes)
}
