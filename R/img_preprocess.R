get_ID <- function (x, label_col='label', return_col='ID'){
        # simplify the series name
        x$series <- gsub ('_cyto_seg.tiff', '', x$series)
        x [, return_col] <- paste (x$series, x [, label_col], sep='_')
        return (x)
}

#' Preprocess the dataset generated from the quantification pipeline
#'
#' @param DAPI a dataframe with information for each nucleus in the DAPI
#' channel
#' @param bright a dataframe with information for each organoid
#' @param tree a dataframe describing the relationship between amnioid and
#' nuclei, required columns: child, parent, condition
#' @importFrom magrittr %>%
#' @export
preprocess <- function (DAPI, bright, tree){
        # Generate a unique label for each object in dataframe
        DAPI <- get_ID (DAPI)
        bright <- get_ID (bright)
        tree <- get_ID (tree, 'child', 'child_ID')
        tree <- get_ID (tree, 'parent', 'parent_ID')

        # clean up the condition name a bit
        #tree$condition <- gsub ('data/', '', tree$condition)
        tree$condition <- tree$condition %>% as.character () %>% basename ()
        #tree$condition <- gsub ('/', '', tree$condition)

        max_DAPI <- DAPI
        # append the information about parent, i.e. amnioid and generate unique
        # labels for them that will match those in `bright`
        add_ID <- tree$parent_ID [match (max_DAPI$ID, tree$child_ID)]
        max_DAPI %>% tibble::add_column (parent_ID=add_ID) %>% 
                dplyr::filter(!is.na (parent_ID)) %>% 
                tidyr::unite ('label_level', c(parent_ID, level), 
                              sep='_', remove=F) -> max_DAPI

        # Generate a unique label for object in each condition in each imaged
        # volumn in each plane
        bright %>%  tidyr::unite ('label_level', c(ID, level), 
                                  sep='_', remove=F) -> bright

        # append the information in `bright` onto `DAPI`, which would be
        # centroids, perimeters, areas, etc
        bright [match (max_DAPI$label_level, bright$label_level), ] %>%
                # make the columns in `bright` distinct from `DAPI`
                magrittr::set_colnames (paste (colnames (bright), '.struct', sep='')) %>%
                cbind (data.frame (max_DAPI)) -> all_data

        ori_data <- all_data
        # prevent duplicate columns from happening
        bind_tree <- tree [, !colnames (tree) %in% c('X', 'series', 'Unnamed..0', 
                                                     'parent_ID')]
        ori_data [match (tree$child_ID, ori_data$ID),] %>% dplyr::select (!condition) %>%
                cbind (bind_tree) %>% dplyr::filter (!is.na (centroid.0.struct)) %>%
                dplyr::select (!all_of (c('X.struct', 'X', 'Unnamed..0.struct',
                        'condition.struct', 'Unnamed..0', 'parent_ID',
                        'child_ID', 'orientation', 'series.struct'))) 
}

#' @importFrom magrittr %>%
#' @export
quantile_plot <- function (dataf, label, y_log=F){
        quan <- seq (0, 1, length.out=100)
        sel_df <- dataf [, label]
        all_quant <- lapply (as.list(1:ncol(sel_df)), function (i){
                                     quantile (sel_df[,i], quan, na.rm=T)})
        all_quant <- do.call (cbind, all_quant)
        all_quant %>% data.frame () %>%
                magrittr::set_colnames (colnames (sel_df)) %>%
                tibble::add_column (x=quan) %>%
                tidyr::gather ('feature', 'value', -x) %>%
                ggplot2::ggplot () + 
                ggplot2::geom_line (ggplot2::aes (x=x, y=value)) +
                ggplot2::facet_wrap (~feature, scales='free_y') +
                TBdev::theme_TB ('dotplot', rotation=0) +
                ggplot2::ylab (label) + ggplot2::xlab ('percentile') -> plot_ob
        if (y_log){
                return (plot_ob + ggplot2::scale_y_log10 ())
        }else{return (plot_ob)}
}

#' Plot distribution of features in amnioid
#' 
#' @export
amnioid_hist <- function (plot_df, xaxis, bin_group=20, convert_radian=F){
        plot_df %>% dplyr::group_by (ID.struct) %>% 
                dplyr::summarise (mean_val = mean (!!as.symbol (xaxis))) -> plot_data
        if (convert_radian){plot_data$mean_val <- plot_data$mean_val/pi*180}
        plot_data$mean_val %>% sort () %>% diff () %>% mean() -> bw
        bw <- bw*bin_group
        ggplot2::ggplot (plot_data) +
                ggplot2::geom_histogram (ggplot2::aes (x=mean_val), binwidth=bw, 
                                         fill='gray', color='black')+
                ggplot2::geom_density (ggplot2::aes (x=mean_val, y=bw*..count..), 
                                       color='red')+
                TBdev::theme_TB('dotplot', rotation=0)+
                ggplot2::xlab (xaxis)
}

#' @importFrom magrittr %>%
#' @export
multi_scatter <- function (plot_df, xaxis, yaxes){
        plot_df %>% dplyr::mutate (xval = !!as.symbol (xaxis)) %>%
                dplyr::select (dplyr::all_of (c('xval', yaxes))) %>%
                tidyr::gather ('feature', 'value', -xval) %>%
                ggplot2::ggplot (ggplot2::aes (x=xval, y=value)) +
                ggplot2::geom_point () +
                ggplot2::geom_smooth (method='lm')+
                ggpubr::stat_cor (ggplot2::aes (label=..rr.label..)) +
                ggplot2::facet_wrap (~feature, scales='free_y') +
                TBdev::theme_TB ('dotplot', rotation=0)+
                ggplot2::xlab (xaxis)
}

#' @importFrom magrittr %>%
#' @export
plot_xy_per_path <- function (x_df, one_path, AP=NULL, genes=c('CGB', 'HLAG'),
                              control_name=c('base', 'OKAE'),
                              condition_col='condition', logxy=F){
        AP <- return_aes_param (AP)
        x_df %>% dplyr::filter (!!as.symbol (condition_col) %in% c(one_path, control_name)) -> plot_data
        plot_data %>% dplyr::group_by (condition) %>%
                dplyr::summarise (mean1=mean (!!as.symbol (genes[1])),
                                  mean2=mean (!!as.symbol (genes[2]))) -> text_df

        ggplot2::ggplot (plot_data) +
                ggplot2::geom_point (ggplot2::aes_string (x=genes[1], 
                        y=genes[2], fill='condition'), color='white',
                            shape=21, size=AP$pointsize)+
                ggrepel::geom_text_repel (ggplot2::aes (x=mean1, y=mean2, 
                                label=condition), data=text_df, size=AP$point_fontsize)+
                TBdev::theme_TB('dotplot', rotation=0,
                                feature_vec=plot_data$condition, aes_param=AP, color_fill=T)+
                ggplot2::ggtitle (one_path) -> plot_ob
        if (logxy){
                plot_ob < plot_ob + ggplot2::scale_x_log10 () +ggplot2::scale_y_log10()
        }
        return (plot_ob)
}

#' @importFrom magrittr %>%
#' @export
boxplot_gene <- function (x_df, one_path, AP=NULL, genes=c('CGB', 'HLAG'),
                          control_name=c('base', 'OKAE'),
                          condition_col='condition', logy=F){
        AP <- return_aes_param (AP)
        x_df %>% dplyr::mutate (pathway=!!as.symbol (condition_col)) %>%
                dplyr::filter (pathway %in% c(one_path, control_name)) -> plot_data
        plot_data %>% dplyr::select (dplyr::all_of (c(genes, 'pathway'))) %>%
                tidyr::gather ('gene', 'expr', -pathway) %>%
                ggplot2::ggplot (aes (x=pathway, y=expr)) +
                ggplot2::geom_boxplot ()+
                ggplot2::facet_wrap (~gene) +
                TBdev::theme_TB ('dotplot', rotation=0, aes_param=AP)+
                ggplot2::ggtitle (one_path) -> plot_ob
        if (logy){
                plot_ob <- plot_ob + ggplot2::scale_y_log10()
        }
        return (plot_ob)
}

wilcox_one_group <- function (xx, comp, ref, expr_col, group_col){
        group1 <- xx [xx[, group_col] == comp, expr_col]
        group2 <- xx [xx[, group_col] == ref, expr_col]
        return (stats::wilcox.test (group1, group2)$p.value)
}

compare_average_one_gene <- function (x_df, group.by, gene, ref='control'){
        groups_mean <- x_df %>% 
                dplyr::select (dplyr::all_of (c(gene, group.by))) %>%
                magrittr::set_colnames (c('gene', 'condition')) %>%
                dplyr::group_by (condition) %>%
                dplyr::summarise (mean_expr = mean (gene, na.rm=T)) %>% 
                data.frame ()
        groups_mean$pval <- sapply (groups_mean$condition, function (x){
                        wilcox_one_group (x_df, x, ref, gene, group.by)})
        groups_mean$feature <- gene
        groups_mean$control <- groups_mean [groups_mean$condition==ref, 'mean_expr']
        groups_mean %>% #dplyr::filter (condition != ref) %>%
                dplyr::relocate (control, .after=mean_expr)
}

#' @export
compare_average <- function (x_df, group.by, genes=c('HLAG', 'CGB', 'TFAP2C'), ref='control'){
        all_dfs <- lapply (as.list(genes), function (x){
                                   compare_average_one_gene (x_df, group.by, x, ref)})
        all_dfs <- do.call (rbind, all_dfs)
        all_dfs$p.adjust <- p.adjust (all_dfs$pval)
        return (all_dfs)
}
