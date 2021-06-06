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

# ----------Plotting----------

#' @importFrom magrittr %>%
#' @export
quantile_plot <- function (dataf, label, y_log=F){
        quan <- seq (0, 1, length.out=100)
        sel_df <- dataf [, label, drop=F]
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
        if (logy){plot_ob <- plot_ob + ggplot2::scale_y_log10() }
        return (plot_ob)
}

remove_na <- function (x){x[is.na (x)] <- 0; return (x)}

get_cond_ord <- function (){
        c('OKAE', 'base', 'CHIR', 'EGF', 'FGF', 'FK', 'PD03')
}

#' Boxplot or Violin plot the count of cells expressing a certain marker
#'
#' @param x_df a dataframe
#' @param gene which gene to plot
#' @param group.by the experimental conditions
#' @param rep.by the replication for each condition
#' @param box_plot whether to use boxplot or violin plot
#' @param express_thres the fluorescence value by which a cell is called
#' positive for a gene
#' @param perc whether to show the percentage of positive cells or the number
#' of positive cells
#' @param add_pval which comparison groups in `group.by` to compare with
#' `control_group`. If it is 'all', all pairs will be compared
#' @importFrom magrittr %>%
#' @export
frequency_plot <- function (x_df, gene, group.by='condition', rep.by='series',
                            plot_type='box', express_thres=15, perc=T, add_pval=NULL,
                            control_group='base', AP=NULL, sum_table=F){
        AP <- return_aes_param (AP)
        x_df %>% dplyr::group_by (!!as.symbol(group.by), !!as.symbol (rep.by)) %>%
                dplyr::mutate (gene_pos = ifelse (!!as.symbol (gene) > express_thres, 
                                                  'positive', 'negative')) %>%
                dplyr::count (gene_pos) %>% dplyr::rename (num=n) %>%
                tidyr::spread (gene_pos, num) %>%
                dplyr::mutate_if (is.numeric, remove_na ) %>% data.frame () -> x_proc

        # order the conditions along the x axis
        x_proc [, group.by] <- partial_relevel (x_proc [, group.by], get_cond_ord () )
        if (perc){
                # for percentage of positive cells
                x_proc %>% dplyr::mutate (freq = 100*positive/(positive + negative) ) -> x_proc
                y_lab <- paste (gene, '+ cells (%)', sep='')
        }else{
                # for the absolute number of positive cells
                x_proc %>% dplyr::mutate (freq = positive) -> x_proc
                y_lab <- paste (gene, '+ cells', sep='')
        }
        if (sum_table){
                x_proc %>% dplyr::group_by (condition) %>% 
                        dplyr::summarise_if (is.numeric, mean) %>% data.frame () %>% print ()
        }
        if (plot_type=='error'){
                x_proc %>% group_by (!!as.symbol (group.by)) %>%
                        dplyr::mutate (y_mean = mean (freq)) %>%
                        dplyr::mutate (y_min = y_mean - sd (freq)) %>%
                        dplyr::mutate (y_max = y_mean + sd (freq)) %>%
                        data.frame () -> x_proc
        }

        ggplot2::ggplot (x_proc, ggplot2::aes_string (x=group.by, y='freq') ) -> plot_ob

        if (plot_type=='box'){plot_ob <- plot_ob +ggplot2::geom_boxplot (ggplot2::aes_string (fill=group.by)) 
        }else if (plot_type=='error'){
                plot_ob <- plot_ob + ggplot2::geom_errorbar (aes_string (
                        color=group.by, ymin='y_min', ymax='y_max'))+
                ggplot2::geom_errorbar (aes_string (color=group.by, ymin='y_mean', ymax='y_mean'))+
                ggplot2::geom_jitter (aes_string (color=group.by), size=AP$pointsize/2) 
        }else{plot_ob <- plot_ob +ggplot2::geom_point (ggplot2::aes_string (fill=group.by), 
                                                       shape=21, size=AP$pointsize*2, color='white') }

        # Ironically for R, appending adjusted p value is complicated
        if (!is.null (add_pval)){
                # only keep interested comparisons
                if (add_pval != 'all'){
                        x_proc %>% dplyr::filter (!!as.symbol (group.by) %in% c(add_pval, 
                                                control_group)) -> x_fil
                }else{x_fil <- x_proc}
                x_fil [, group.by] <- partial_relevel (x_fil [, group.by], get_cond_ord () )
                x_fil %>% ggpubr::compare_means (as.formula (paste ( 'freq ~', group.by)), 
                                       data=., ref.group=control_group) -> pstat
                # calculate the position for the p value labels
                rstatix::get_y_position (data=x_fil, formula=as.formula (paste ( 'freq ~', group.by)), 
                                         ref.group=control_group) -> y_pos

                pstat$y.position <- y_pos$y.position
                plot_ob <- plot_ob + ggpubr::stat_pvalue_manual(pstat,
                label='p.adj', size=AP$point_fontsize) 
        }
        plot_ob+theme_TB ('dotplot', color_fill=T, rotation=90, aes_param=AP)+
                #ggplot2::guides (fill=ggplot2::guide_legend ())+ 
                ggplot2::ylab (y_lab)+ggplot2::theme (legend.position='none')
}

fluoro_plot <- function (x_df, gene, group.by='condition',
                         control_group='base', color.by=NULL, color_thres=15,
                         AP=NULL){
        AP <- return_aes_param (AP)
        AP$color_vec <- c(AP$color_vec, 'unknown_gray'='#bebebe')
        x_df [, group.by] <- partial_relevel (x_df [, group.by], get_cond_ord () )
        x_df %>% group_by (!!as.symbol (group.by)) %>%
                dplyr::summarise(y_mean = mean (!!as.symbol (gene)), 
                                 y_sd = mean (!!as.symbol (gene))) %>%
                dplyr::mutate (y_min = y_mean - y_sd, y_max= y_mean + y_sd) %>%
                data.frame () -> x_sum

        if (is.null (color.by) & !is.null (color_thres)){
                cell_type <- c(EVT='HLAG', STB='CGB', CTB='TFAP2C')
                tissue_color <- names (cell_type) [cell_type==gene]
                color.by <- 'Type'
                x_df %>% dplyr::mutate (Type = ifelse (!!as.symbol (gene)> color_thres, 
                                                       tissue_color, 'unknown')) -> x_df
        }

        x_df %>% ggpubr::compare_means (as.formula (paste (gene, '~', group.by)), 
                               data=., ref.group=control_group) -> pstat
        pstat$y.position <- max (x_df [, gene])

        ggplot2::ggplot (x_df, aes_string (x=group.by))+
                ggplot2::geom_jitter (aes_string (y=gene, color=color.by), 
                                      size=AP$pointsize/4, height=0) +
                ggplot2::geom_errorbar (aes_string (
                        ymin='y_min', ymax='y_max'), data=x_sum)+
                ggplot2::geom_errorbar (aes_string (ymin='y_mean', ymax='y_mean'), data=x_sum)+
                ggpubr::stat_pvalue_manual(pstat, label='p.adj', x='group2',
                                           size=AP$point_fontsize*0.5) +
                theme_TB ('dotplot', feature_vec=x_df [, color.by], aes_param=AP)+
                ggplot2::theme (legend.position='none')
}

# ----------Statistics----------

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

# ----------clonogenicity----------

#' Count the number of colonies
#'
#' @description This function first applies the following filter: remove DBSCAN
#' outliers, area > 400, number of cells per colony outside the range of 6~300
#' Secondly, it counts the number of colonies that meet such criteria.
#' @param dat a dataframe containing the following columns
#' `series`: area ID in a culture well
#' `Date`: image batch
#' `cluster`: the cluster ID assigned to each `series`
#' `condition`: each siRNA condition
#' @return a dataframe with the following columns:
#' `series`, `Date`, `condition`
#' `n`: number of colonies
#' `rel_n`: colony number relative to the mean in the control condition
#' @importFrom magrittr %>%
#' @export
colony_count <- function (dat, control_condition='GFP', cell_lab=NULL){
        dat %>% tidyr::unite ('series', c('series', 'Date'), sep='__') %>%
        dplyr::filter (cluster != -1) %>% 
        dplyr::filter (area < 400) %>%
        dplyr::count (condition, series, cluster) %>%
        dplyr::filter (n>=6 & n < 300) %>% 
        dplyr::rename (num_cell= n) %>%
        dplyr::group_by (condition, series) %>%
        dplyr::summarise (num_colony= dplyr::n(), mean_colony=mean (num_cell) ) %>%
        data.frame () %>%
        tidyr::separate (series, c('series', 'Date'), sep='__') -> proc_dat
        if (!is.null (control_condition)){
                all_dates <- unique (proc_dat$Date)# all siRNA batches
                proc_dat$rel_num_colony <- proc_dat$num_colony
                proc_dat %>% dplyr::filter (condition == control_condition) %>%
                        dplyr::group_by (Date) %>%
                        dplyr::summarise (mean_con_n = sum(num_colony)) %>%
                        data.frame() -> control_dat 
                print (control_dat)
                for (i in all_dates){
                        one_cond <- proc_dat$Date == i
                        mean_control_n <- control_dat [control_dat$Date==i, 'mean_con_n']
                        proc_dat$rel_num_colony [one_cond] <- proc_dat$rel_num_colony [one_cond]/mean_control_n
                }
        }
        if (is.null (cell_lab)) {data (siRNA_cell, package='TBdev'); cell_lab <- siRNA_cell}
        proc_dat$celltype <- cell_lab$celltype [match (proc_dat$condition, toupper (cell_lab$TF))]
        return (proc_dat)
}

#' Boxplot of the colony features in each siRNA condition
#'
#' @param dat a dataframe from `colony_count`
#' @param y_val y axis of the boxplot
#' @param control_condition which is the control group in the `condition`
#' column
#' @param plot_pval which pvalue format to plot, can be 'p.adj' for the BH
#' adjusted p value score, 'p.signif' for the symbol for 'p.adj' or 'p' for
#' unadjusted p value
#' @export
plot_colony_count <- function (dat, y_val='rel_num_colony', shape_by='Date',
                               control_condition='GFP', plot_pval='p.signif',
                               plot_type='box', paired_test=F, display_shape=T,
                               new_level=NULL, AP=NULL){
        AP <- return_aes_param (AP)
        AP$color_vec <- c('control'='#32CD32', AP$color_vec)
        dat %>% dplyr::select(celltype, condition) %>%
                dplyr::filter (!duplicated (condition)) %>%
                tibble::deframe () -> cond_cell
        ord_cond <- partial_relevel (names (cond_cell))
        cond_cell <- cond_cell [order(ord_cond)]
        dat$condition <- factor (dat$condition, levels=as.character (cond_cell))
        if (!is.null (new_level)){dat$condition <- partial_relevel (dat$condition, new_level)}

        dat %>% dplyr::filter (condition==control_condition) %>% 
                dplyr::summarise (mean_n = median(!!as.symbol (y_val))) %>% 
                tibble::deframe () -> mean_lev

        if (!paired_test){
                ggpubr::compare_means (as.formula (paste (y_val, '~condition', sep='')), 
                                       data=dat, ref.group=control_condition,
                                       p.adjust.method='BH') -> stat_test
        }else{
                pairwise_p (dat, control_condition=control_condition,
                            y_val=y_val, group.by=shape_by) -> stat_test
        }
        if (shape_by %in% colnames (dat)){
                dat %>% data.frame () -> dat
                dat [, shape_by] <- as.character (dat [, shape_by])
        }
        max_y <- max (dat [, y_val])
        if (plot_type=='error'){
                dat %>% group_by (condition) %>%
                        dplyr::mutate (y_mean = mean (!!as.symbol (y_val))) %>%
                        dplyr::mutate (y_min = y_mean - sd (!!as.symbol (y_val))) %>%
                        dplyr::mutate (y_max = y_mean + sd (!!as.symbol (y_val))) %>%
                        data.frame () -> dat
        }
        ggplot2::ggplot (dat, aes_string (x='condition', y=y_val))+
                ggplot2::geom_hline (yintercept=mean_lev, color='red', linetype='dashed') +
                ggpubr::stat_pvalue_manual (stat_test, y.position=max_y,
                                            label=plot_pval, xmin='group2',
                                            xmax='group2', tip.length=0,
                                            label.size=AP$point_fontsize,
                                            angle=30)+
                TBdev::theme_TB ('dotplot', rotation=45, feature_vec=dat$celltype, aes_param=AP)+
                ggplot2::theme (aspect.ratio=0.6) -> plot_dat
                
        if (plot_type=='box'){plot_dat+
                ggplot2::geom_boxplot (aes (color=celltype), show.legend=F) %>% return ()
        }else if (plot_type=='error') {
                plot_dat+ ggplot2::geom_errorbar (aes (
                        color=celltype, ymin=y_min, ymax=y_max))+
                ggplot2::geom_errorbar (aes (color=celltype, ymin=y_mean, ymax=y_mean))+
                ggplot2::geom_jitter (aes_string (color='celltype'), 
                                      size=AP$pointsize/2) %>% return ()
        }else{
                if (display_shape){
                        return (plot_dat+ggplot2::geom_point (aes_string (color='celltype', 
                                        shape=shape_by), size=AP$pointsize*2))
                }else{
                        return (plot_dat+ggplot2::geom_point (aes_string (color='celltype', 
                                        ), size=AP$pointsize*2))
                }
}}

pairwise_p <- function (dat, control_condition='GFP', y_val='total_num', group.by='Date'){
        dat %>% data.frame () -> dat
        dat %>% dplyr::pull (condition) %>% unique () %>% as.list () %>%
        lapply (function (xx){
                dat %>% dplyr::filter (condition %in% control_condition) -> control 
                dat %>% dplyr::filter (condition %in% xx) -> experiment
                control <- control [match (experiment[, group.by], control[, group.by]), ]
                stats::t.test (control [, y_val], experiment [, y_val], paired=T) -> test_stat
                data.frame('.y.'=y_val, group1=control_condition, group2=xx, p=test_stat$p.value)}) %>% 
        do.call(what=rbind) -> stat_dat
        stat_dat$p.adj <- stats::p.adjust (stat_dat$p, method='BH')
        stat_dat$p.signif <- stats::symnum(as.numeric (stat_dat$p.adj),
                cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1), na='',
                symbols = c("***", "**", "*", ".", "ns"))
        return (stat_dat)
}

#' @export
reproducibility_mat <- function (dat, count_col='total_num'){
        all_levels <- unique (dat$Date)
        N <- length(all_levels)
        comp_mat <- matrix (0, N, N)
        for (i in 1:N){
                for (j in 1:N){
                        if (i != j){
                                i_dat <- dat [dat$Date==all_levels[i],]
                                j_dat <- dat [dat$Date==all_levels[j],]
                                if (nrow (i_dat) < nrow (j_dat)){
                                        j_dat <- j_dat [match (i_dat$condition, j_dat$condition), ]
                                }else{
                                        i_dat <- i_dat [match (j_dat$condition, i_dat$condition), ]
                                }
                                comp_mat [i, j] <- stats::cor (i_dat [, count_col], 
                                                               j_dat [, count_col])
                        }else{
                                comp_mat [i, j] <- 1
                        }
                }
        }
        colnames (comp_mat) <- all_levels
        rownames (comp_mat) <- all_levels
        return (comp_mat)
}
