# This script processes the raw data generated from BRGP algorithm.
# workflow: seurat: get expression matrix > gpflow: GPLVM > STREAM: graph based
# trajectory inference > gpflow: BRGP > tidyverse: statistics

# ----------Plotting----------

#' Plot pseudotimee against real time
#'
#' @param dat a dataframe
#' @param real_time the column in `dat` that stores the real time information
#' @param pseudotime the column in `dat` that stores the pseudotime information
#' @param color_by the column that stores the features to colors
#' @importFrom ggplot2 aes_string
#' @importFrom magrittr %>%
#' @export
pseudo_real_time <- function (dat, real_time, pseudotime, color_by, AP=NULL){
        AP <- return_aes_param (AP)
        gsub ('^D', '', dat[, real_time]) %>% as.numeric () -> numer_date
        time_cor <- stats::cor (numer_date, dat [, pseudotime])
        dat %>% ggplot2::ggplot (aes_string (x=real_time, y=pseudotime, fill=color_by) ) + 
                ggplot2::geom_point (shape=AP$normal_shape, color=AP$point_edge_color, 
                                     size=AP$pointsize)+
                ggplot2::ggtitle ( paste ('\u03c1 =', format (round (time_cor, 2), nsmall=2) ) ) +
                theme_TB ('dotplot', feature_vec=dat [, color_by], color_fill=T, rotation=90, AP=AP)+
                custom_tick (dat [, pseudotime]) +
                ggplot2::labs (fill='') + ggplot2::ylab ('pseudotime') 
}

#' Remove the `mean_` and `var_` prefixes from a vector
get_gene_names <- function (x){
        x_no_mean <- gsub ('^mean_', '', x)
        x_no_var  <- gsub ('^var_' , '', x_no_mean)
        return (x_no_var)
}

#' Return parts of the dataframe `x_df` that begins with `sel_col` and with
#' rows being `genes`
preprocess_df <- function (x_df, sel_col, genes){
        sel_column <- paste (sel_col, genes, sep='')
        sel_column <- sel_column [ sel_column %in% colnames (x_df) ]
        x_df <- x_df [, c(sel_column, 'x', 'branch')]
        colnames (x_df) <- get_gene_names (colnames (x_df) )
        x_df$info_type <- sel_col
        return (x_df )
}

#' Gene expression with uncertainty over pseudotime
#' 
#' @param x dataframe containing mean + var information estimated by BRGP
#' @param exp_mat the expression matrix, i.e. true expression levels
#' @param genes the genes to be plotted
#' @importFrom ggplot2 aes aes_string
#' @importFrom magrittr %>%
#' @export
gene_over_pseudotime <- function (x, exp_mat, genes, metadata, color_feature,
                                  num_col=4, num_row=NULL, branch_assignment=NULL, 
                                  peak_data=NULL, gene_col='feature', AP=NULL){
        AP <- return_aes_param (AP)
        mean_df <- preprocess_df (x, 'mean_', genes)
        var_df  <- preprocess_df (x, 'var_' , genes)
        join_df <- rbind (mean_df, var_df)
        print ('getting data for ribbon plot')
        new_genes <- as.character (genes [ genes %in% colnames (exp_mat) ])
        join_df %>% tidyr::gather ('gene', 'val', -x, -branch, -info_type) %>%
                tidyr::spread (info_type, val) %>% 
                dplyr::mutate (ymin = mean_ - 2*sqrt (var_) ) %>%
                dplyr::mutate (ymax = mean_ + 2*sqrt (var_) ) %>%
                dplyr::mutate (gene = factor (gene, levels=new_genes) )-> plot_df

        if (!is.null(branch_assignment)){
                print ('reassigning branch names')
                plot_df$branch <- branch_assignment [as.factor (plot_df$branch) ]
        }
        print ('processing raw expression matrix')
        exp_mat %>% tibble::add_column (cell_ID = rownames (exp_mat) ) %>%
                dplyr::select (c(new_genes, 'pseudotime', 'cell_ID')) %>%
                tidyr::gather ('gene', 'mean_', -pseudotime, -cell_ID) %>%
                dplyr::mutate (gene = factor (gene, levels=new_genes) ) -> point_df
        point_df$color_by <- metadata [match (point_df$cell_ID, rownames (metadata)), color_feature]
        point_df %>% dplyr::filter (!is.na (color_by)) -> point_df

        print ('plotting')
        ggplot2::ggplot (plot_df ) +
                ggplot2::geom_point (aes (x=pseudotime, y=mean_, color=color_by), data=point_df, shape=20)+
                ggplot2::geom_ribbon (aes (x=x, y=mean_, ymin=ymin, ymax=ymax, fill=branch), alpha=0.8 )+
                ggplot2::facet_wrap (~gene, scales='free', ncol=num_col, nrow=num_row) +
                theme_TB ('no_arrow', feature_vec=point_df$color_by, AP=AP) +
                theme_TB (feature_vec = plot_df$branch, AP=AP, color_fill=T) +
                ggplot2::xlab ('pseudotime') + ggplot2::ylab ('gene expression') +
                ggplot2::labs (color='cell type') -> plot_ob

        print ('plot vertical lines at peak times')
        if (!is.null(peak_data)){
                selected_peak <- peak_data [peak_data [, gene_col] %in% genes,]
                plot_ob <- plot_ob + plot_peak_line (selected_peak, AP, gene_col)
        }
        return (plot_ob)
}

#' Heatmap shows gene expression along the order of pseudotime across groups
#' 
#' @description I did not find this function particularly helpful.
#' `seurat_heat` is more flexible.
pseudotime_heat <- function (seurat_list, show_genes, group.by, order_row='pseudotime',
                             assay='RNA', slot_data='counts', return_sep=F, col_title=NULL,
                             ...){
        heat_list <- list()
        N <- length (seurat_list)
        if (is.null (col_title) ){col_title <- names (seurat_list) }
        if (is.null (col_title) ){col_title <- rep (NULL, N)}
        for (i in 1:N){
                one_seurat <- seurat_list[[i]]
                assay_dat <- Seurat::GetAssayData (one_seurat, assay=assay, slot=slot_data)
                PT_order <- order ( assay_dat [order_row,])
                if (i==N){show_col_anna <- T}else{show_col_anna <- F}
                heat_list [[i]] <- seurat_heat (one_seurat, color_row=show_genes, 
                                                    group.by = group.by,
                                 column_rotation=0, 
                                 slot_data =slot_data, assay=assay, 
                                 column_title=col_title [i],
                                 group_order=PT_order, 
                                 show_column_anna = show_col_anna,
                                 column_split=NULL, ...)
        }
        if (return_sep){return (heat_list)
        }else{
                final_plot <- heat_list[[1]]
                for (i in 2:N){
                        final_plot <- final_plot + heat_list[[i]]
                }
                return (final_plot)
        }

}

#' Convert matrix into a list of seurat objects based on branch assignment
#'
#' @param exp_mat a matrix with rows being the cells and columns being the
#' genes. There must be rownames and colnames
#' @param ref_meta a dataframe containing information for the cells in `exp_mat`
#' @param branch_lab a vector with the same length as the row number of
#' `exp_mat`, containing the branch assignment of each cell
#' @param label_list per item of returned seurat object in the return value
#' list, cells with which labels are chosen. It should be a list
#' @return a list of seurat objects of different branch assignments
mat_to_seurat <- function (exp_mat, ref_meta, branch_lab, label_list){
        print ('assigning metadata including branch information')
        if (!is.null(ref_meta)){
                PT_seurat_meta <- ref_meta [ match (rownames (exp_mat), rownames (ref_meta) ), ] 
        }else{ 
                PT_seurat_meta <- data.frame (ID=1:nrow (exp_mat) )
                rownames (PT_seurat_meta) <- rownames (exp_mat)
        }
        PT_seurat_meta$branch <- branch_lab [ match (rownames (PT_seurat_meta), 
                                                     names (branch_lab)) ]
        print ('converting to seurat object')
        PT_seurat <- Seurat::CreateSeuratObject (t(exp_mat), meta.data=PT_seurat_meta)
        seurat_list <- list ()
        for (i in 1:length(label_list)){
                print (paste ('getting the', i, 'item') )
                seurat_list [[i]] <- PT_seurat [, branch_lab %in% label_list[[i]] ]
        }
        return (seurat_list)
}

#' Load pseudotime information to Seurat
#'
#' @description Obtain the mean values and store them in Seurat objects
#' according to the branch assignment indicated by `label_list`
#' @param pred_all dataframe generated from pseudotime analysis. It contains
#' the mean and variance of all the genes, starting with 'mean_' and
#' 'variance_' respectively. In addition, it has a column called `x` containing
#' pseudotime information and a column called `branch` for branch assignment
#' @param label_list per item of returned seurat object in the return value
#' list, cells with which labels are chosen. It should be a list
#' @return a list of seurat objects of different branch assignments
#' @importFrom magrittr %>%
#' @export
raw_to_seurat <- function (pred_all, label_list){
        pred_list <- sep_mean_val (pred_all)
        pred_list [[1]] %>% data.frame () %>% tibble::add_column (pseudotime=pred_all['x']) %>% 
                as.matrix () -> exp_mat_inf 
        pred_all %>% dplyr::select (branch) %>% tibble::deframe () -> meta_data

        rownames (exp_mat_inf) <- paste ('cell', 1:nrow(exp_mat_inf), sep='')
        names (meta_data) <- rownames (exp_mat_inf)
        return (mat_to_seurat (exp_mat_inf, NULL, meta_data, label_list))
}

#' Integrate expression matrix, branch information, pseudotime data 
#'
#' @description I did not find this function helpful because later I added them
#' to the metadata of the integrated dataset.
integrated_seurat <- function (exp_mat, pseudotime, branch, metadata){
        data (all_types, package='TBdev')
        int_meta <- data.frame ('pseudotime'=pseudotime, 'branch'=branch)
        assign_branch <- c('main', 'EVT_branch', 'STB_branch')
        int_meta$branch <- assign_branch [ as.factor (int_meta$branch) ]
        rownames (int_meta) <- colnames (exp_mat)
        int_meta <- cbind ( int_meta, metadata [ match ( rownames (int_meta), rownames (metadata) ), ]  )
        inter_seurat <- Seurat::CreateSeuratObject (exp_mat, meta.data=int_meta )
        inter_seurat <- inter_seurat [, !(inter_seurat$broad_type %in% all_types$non_TB_lineage) ]
        return (inter_seurat)
}

#' Calculate the pathway module score for the Seurat object for each branch
#'
#' @param seurat_list a list of seurat objects
#' @param save_dir where the module scores are saved
seurat_list_score <- function (seurat_list, save_dir, label='DF_B', KeggID=NULL){
        module_list <- list ()
        if (is.null (KeggID)){data (KeggID, package='TBdev')}
        for ( i in 1:length (seurat_list) ){
                save_path <- paste (save_dir, paste (label, i, '.csv', sep=''), sep='/')
                module_score <- get_module_score (seurat_list[[i]], all_path=KeggID, 
                                                  save_path= save_path)
                colnames (module_score) <- rownames (seurat_list[[i]]@meta.data)
                module_score ['pseudotime', ] <- as.vector (Seurat::GetAssayData (
                                seurat_list[[i]], slot='counts') ['pseudotime', ])
                module_list [[i]] <- Seurat::CreateSeuratObject ( module_score,  
                                                meta.data=seurat_list[[i]]@meta.data)
        }
        return (module_list)
}

# ----------PCA----------

#' PCA plot of cells with pseudotime trajectory projected onto it
#'
#' @param pc_pt a prcomp object
#' @param pt_mat the matrix for pseudotime trajectory
#' @param branch_info a vector for branch assignment
#' @description The reason why I have not integrated PCA computation into this
#' function is that pca takes too long for large matrices. That's why the
#' function looks quite inconvenient. I use `gmodels::fast.prcomp` for PCA.
#' Similarly, users need to calculate the PC projected variance externally e.g.
#' proj_var <- solve (pc_pt$rotation^2, t(pred_list[[2]]), tol=1e-20)
#'
#' Strictly speaking, projecting variance onto PCs require NNLS. I have
#' implemented it in python but no time to do so in R. Fortunately, I did not
#' need to use this function often.
#' @importFrom ggplot2 aes_string
#' @export
pca_with_pt_line <- function (pc_pt, pt_mat, metadata, color.by,
                              branch_info=NULL, proj_var=NULL, num_dim=c(1,2),
                              AP=AP){
        AP <- return_aes_param (AP)
        rotated <- pt_mat %*% pc_pt$rotation
        plot_data <- data.frame (pc_pt$x [, num_dim])
        plot_data <- cbind (plot_data, metadata [match (rownames (plot_data), rownames (metadata) ), ] )
        if (is.null (branch_info)){branch_info <- rep ('branch0', nrow(rotated) ) }
        plot_line <- data.frame (rotated [, num_dim]) %>% tibble::add_column (branch=branch_info)

        x_axis <- colnames (plot_data)[1]
        y_axis <- colnames (plot_data)[2]
        ggplot2::ggplot (plot_data, aes_string (x=x_axis, y=y_axis) ) + 
                ggplot2::geom_point (aes_string (fill= color.by), shape=AP$normal_shape, 
                                     size=AP$pointsize, color=AP$point_edge_color) +
                ggplot2::geom_line (aes_string (x=x_axis, y=y_axis, color='branch'), data=plot_line, size=2) -> plot_ob
        if (!is.null (proj_var)){
                var_data <- pmax (t(proj_var) [, num_dim], 0) %>% data.frame () 
                min_val <- plot_line [, y_axis] - 2*sqrt (var_data [, y_axis])
                max_val <- plot_line [, y_axis] + 2*sqrt (var_data [, y_axis])
                var_data [, 'branch'] <- plot_line [, 'branch']
                var_data [, 'ymin'] <- min_val 
                var_data [, 'ymax'] <- max_val 
                var_data [, 'x'] <- plot_line [, x_axis]
                plot_ob <- plot_ob + ggplot2::geom_ribbon (aes_string (x='x', ymin='ymin', 
                                                              ymax='ymax'), fill='gray', data=var_data)
        }
        plot_ob + theme_TB('dim_red', plot_ob=plot_ob, feature_vec = 
                           plot_data [, color.by], color_fill=T, AP=AP)
}


#' Calculate mean log likelihood
#'
#' @description This function works extremely slow. I recommend using the
#' python alternative.
mean_log_likelihood <- function (mu, variance, x){
    out = -(mu^2 + colMeans (x^2) -2*mu*colMeans (x))/variance*2 
    out = out - sqrt(variance) - 0.5*log (2*pi)
    return (exp (rowMeans (out) ))
}

# ----------temporal clustering----------

#' Plot the clustering results of times against the peak gradient
#'
#' @param peak_plot a dataframe generated by `find_peak` or
#' `find_transition_from_seurat`. It must contain the column `peak_time` and `val`,
#' the latter for selection of the most significant genes
#' @param metadata a dataframe for adding the cell types colors
#' @param color_by which column in `peak_plot` to color the points
#' @param color_bar which column in `metadata` contains cell type information
#' to provide a cellular context for the pseudotime line
#' @param time_col which column in `metadata` contains pseudotime information
#' for each cell
#' @param exclude_perc by how much percent the pseudotime line should be
#' truncated. This is because pseudotime points sometimes contain some outliers
#' that may create a lot of white space.
#' @importFrom ggplot2 aes aes_string
#' @importFrom stats quantile
#' @importFrom magrittr %>%
#' @export
time_cluster_plot <- function (peak_plot, metadata, show_text_prop=0.95,
                               color_by='cluster', color_bar='broad_type',
                               time_col='MGP_PT', exclude_perc=0.02, vjust=0.,
                               thickness_ratio=0.05, repel_force=1,
                               repel_point=NULL, AP=NULL){
        AP <- return_aes_param (AP)
        PT_type <- data.frame (PT=metadata [, time_col], Type=metadata [, color_bar]) 
        min_val <- quantile(peak_plot$peak_time, exclude_perc)
        max_val <- quantile(peak_plot$peak_time, 1-exclude_perc)
        PT_type %>% dplyr::filter (PT > min_val & PT < max_val) -> PT_type

        peak_plot %>% dplyr::group_by (!!as.symbol (color_by) ) %>%
                dplyr::filter (val > quantile (val, show_text_prop) ) %>% 
                dplyr::arrange (dplyr::desc(val)) -> label_peak

        if (!is.null (repel_point)){
                peak_plot %>% dplyr::group_by (!!as.symbol (color_by) ) %>%
                        dplyr::filter (val > quantile (val, repel_point) ) %>% 
                        dplyr::arrange (dplyr::desc(val)) -> new_peak
                new_peak <- new_peak [!new_peak$feature %in% label_peak$feature,]
                new_peak [, color_by] <- NA
                rbind (label_peak %>% dplyr::ungroup(), 
                       new_peak %>% dplyr::ungroup()) -> label_peak
        }

        min_y <- min (peak_plot [, 'val']) - vjust
        max_y <- max (peak_plot [, 'val']) - vjust
        thickness <- (max_y - min_y)*thickness_ratio

        ggplot2::ggplot (peak_plot, aes (x=peak_time, y=val) ) + 
                ggplot2::geom_point (aes_string (color=color_by), size=AP$pointsize) + 
                ggrepel::geom_text_repel (aes_string (label='feature', color=color_by), 
                                          data=label_peak, show.legend=F, force=repel_force, 
                                          size=AP$point_fontsize, fontface='bold') + 
                theme_TB ('dotplot', feature_vec = as.character (peak_plot [, color_by]), 
                          color_fill=F, rotation=0, AP=AP)+
                ggplot2::geom_ribbon (aes_string (x='PT', fill='Type', ymax=min_y, ymin=min_y-thickness),
                            size=AP$pointsize, data=PT_type, inherit.aes=F) +
                add_custom_color (feature_vec = PT_type$Type, aes_param=AP, color_fill=T)+
                custom_tick (peak_plot$val) +
                ggplot2::theme (aspect.ratio=0.5) + ggplot2::labs (fill =color_bar) + 
                ggplot2::ylab ('maximum gradient') + ggplot2::xlab ('pseudotime') +
                ggplot2::xlim ( c(min_val, max_val) )
}

