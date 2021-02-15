# devtools::install_github('immunogenomics/presto')

#' filter out genes with low expressions
#' 
#' @description I am not really sure the most efficient way to do this given
#' that Seurat object does not have rowData. Therefore this function is quite
#' slow.
#' @param x a Seurat object
#' @param quantile_val filter away genes whose expressions are below a certain
#' percentile
#' @export
filter_genes <- function (x, quantile_val, assay='RNA', slot_data='data'){
        exp_mat <- as.matrix (Seurat::GetAssayData (x, assay=assay, slot=slot_data) )
        mean_gene <- rowMeans (exp_mat)
        filter_thres <- stats::quantile (mean_gene, quantile_val)
        filter_genes <- mean_gene > filter_thres
        print (paste ('filtering away', sum (!filter_genes), 'genes'))
        return (x [filter_genes, ])
}

# ----------------
# DE gene analysis
# ----------------

#' Find DE genes across a particular feature
#' 
#' @param x a Seurat object
#' @param group.by a column in the meta data to group the cell types
#' @param directory where to save the computed markers
#' @param filename the name of the data file. If NULL, it will be generated
#' automatically in the format of '$directory/marker_list_$group.by_$label.csv'
#' @param method either 'group' to run DE gene of one group with respect to the
#' rest of the dataset using `presto::wilcoxauc`, or 'pairwise' to run DE gene
#' between any possible pairs of groups using `find_DE_pairwise`. The pairwise
#' option can be very computationally expensive.
#' @importFrom Seurat Idents<-
#' @export
find_DE_genes <- function (x, directory=NULL, group.by='seurat_clusters',
                           filename= NULL, label='all', method='group',
                           assay='RNA', slot_data='data'){

        if (is.null (filename)){
                filename <- paste ('marker_cluster_', group.by, '_', label, '.csv', sep='')
        }
        if (is.null (directory)){proceed <- T
        }else{
                file_name <- paste (directory, filename, sep='/')
                if (!file.exists (file_name)){proceed <- T
                }else{proceed <- F}
        }
        if (proceed){
                Idents (x) <- x@meta.data [, group.by]
                # Don't use the `FindMarkers` function in Seurat which is too
                # slow
                if (method == 'group'){
                        markers <- presto::wilcoxauc (x, group.by, 
                                assay= slot_data, seurat_assay=assay)
                }else if (method == 'pairwise'){
                        markers <- find_DE_pairwise (x, group.by, 
                                slot_data= slot_data, assay=assay)
                }
                if (!is.null (directory)) {utils::write.csv (markers, file_name)}
        }else{
                print ('loading precomputed DE genes')
                # use data table to read csv is much faster for large files
                markers <- data.frame (data.table::fread(file_name))
        }
        return (markers)
}

find_DE_onepair <- function (x, group.by, group1, group2, assay='RNA', slot_data='data'){
        sel_x <- x [, x@meta.data [, group.by] %in% c(group1, group2)]
        onepair_marker <- presto::wilcoxauc (sel_x, group.by, assay=slot_data, 
                                             seurat_assay=assay)
        compare_vec <- onepair_marker$group == group1
        onepair_marker$compare_group <- data.table::fifelse (compare_vec, group2, group1)
        return (onepair_marker)
}

#' Find pairwise DE gene expression
#'
#' @param x a Seurat object
#' @param group.by a column in the meta data to group the cell types
#' @importFrom magrittr %>%
find_DE_pairwise <- function (x, group.by, assay='RNA', slot_data='data'){
        all_types <- x@meta.data [, group.by] %>% unique () %>% as.character ()
        N <- length (all_types)
        all_pairs <- matrix (0, N, N)
        # do not select the diagonal terms
        all_pairs [lower.tri (all_pairs, diag=F)] <- 1
        all_list <- list ()
        k <- 1
        for (i in 1:N){
                for (j in 1:N){
                        if (all_pairs [i, j] == 1){
                                print (paste ('comparing', all_types[i], all_types[j] ))
                                all_list [[k]] <- find_DE_onepair (x, group.by, 
                                all_types[i], all_types[j], assay=assay,
                                slot_data=slot_data)
                                k <- k+ 1
                        }
                }
        }
        return (do.call (rbind, all_list))
}

# ---------------
# Select DE genes
# ---------------

#' Select the top DE genes from each cluster
#'
#' @description If 1 DE gene appears in 2 clusters, the cluster with the
#' highest avg_logFC will be chosen. 
#' Then the user can specify the option in `iterative` whether the other DE
#' genes from the group will be chosen to ensure at the end the number of DE
#' genes in each group is the same.
#' NB: only upregulated genes will be selected
#' 
#' @param DE_genes a dataframe
#' @param top_number how many DE genes per cluster/group to choose
#' @param gene which column has the gene symbols
#' @param weighting which column has the gene rankings
#' @param cluster which column has the cluster names/numbers
#' @param iterative how many iterations to run until the number of DE genes in
#' each group is the same. By default, no iteration will be performed. I
#' recommend choosing a moderate size e.g. 10~20 to ensure the end result. I
#' could have implemented a while loop but until I consider all possible
#' scenarios, I would not do so. 
#' @importFrom magrittr %>%
#' @export
unique_DE_genes <- function (DE_genes, top_number, gene='feature',
                             weighting='logFC', cluster='group', iterative=NULL){
        if (is.null(iterative)){iterative<-1}
        top_num <- top_number
        for (i in 1:iterative){
                DE_frame <- unique_DE_genes_one_iter (DE_genes, top_num,
                                                      gene, weighting, cluster)
                DE_frame %>% dplyr::group_by (!!as.symbol (cluster)) %>%
                        dplyr::slice_max (!!as.symbol(weighting), n=top_number) -> DE_frame
                count_mean <- DE_frame %>% dplyr::count (!!as.symbol (cluster)) %>%
                        tibble::deframe () %>% mean ()
                if (count_mean != top_number){
                        top_num <- top_num + 1
                }else{break}
        }
        return (DE_frame)
}

#' @importFrom magrittr %>%
unique_DE_genes_one_iter <- function (DE_genes, top_number, gene, weighting,
                                      cluster){
        DE_genes %>% dplyr::group_by (!!as.symbol (cluster) ) %>%
                dplyr::top_n (n=top_number, wt=!!as.symbol (weighting) ) %>%
                dplyr::ungroup () %>%
                dplyr::arrange (desc (!!as.symbol (weighting) ) ) %>%
                dplyr::distinct (!!as.symbol (gene), .keep_all=T ) 
}

list_DE_genes <- function (x, directory, group.by, top_number=8, label='all',
                           gene='feature', cluster='group'){
        DE_genes <- find_DE_genes (x, directory, group.by, label=label)
        top_markers <- as.data.frame (unique_DE_genes (DE_genes, top_number))
        DE_genes <- as.vector (top_markers [, gene])
        names (DE_genes) <- top_markers[,cluster]
        return (DE_genes)
}

#' Merge DE genes and known lineage markers
#'
#' @export
DE_lineage_genes <- function (x, directory, top_number=3, label='all',
                              group.by='revised', lineage_markers=NULL){
        if (is.null (lineage_markers)){data (lineage_markers, package='TBdev') }
        DE_genes <- list_DE_genes (x, directory, group.by=group.by,
                                   top_number=top_number, label=label)
        all_genes <- c(lineage_markers,  DE_genes)
        all_genes <- all_genes [ gtools::mixedorder (names (all_genes)) ]
        return (all_genes)
}


# ----------------------------------
# Exporting DE genes for publication
# ----------------------------------

#' Save Differentially Expressed Genes
#'
#' @description This function intends to save DE genes for presentation and
#' publication but not for computing. For downstream analysis, the
#' `find_DE_genes` function would be able to save the raw output.
#' @param x a Seurat object
#' @param save_dir where to save the results
#' @param group.by a column in the meta data to group the cell types
#' @param label the file name would be: '$save_dir/DE_$group.by_$label.$suffix'
#' @param save_format can be 'excel' or 'gene_table' for average expression only
#' @param weighting which column to arrange the genes from top to bottom
#' @param gene which column stores the gene names
#' @param value which column stores the average expression
#' @param show_num how many DE genes to save per group. If 'all', everything
#' will be saved
#' @param AP a list that contains the `cell_order` entry to arrange the cells
#' @param organism_db a genome database, e.g. org.Hs.eg.db. It is optional
#' unless you would like to add a ensemble ID column
#' @seealso `find_DE_genes`
#' @importFrom magrittr %>%
#' @export
save_DE_genes <- function (x, save_dir, markers=NULL, 
                           # save path related
                           group.by='revised', 
                           label='all', 
                           save_format='excel',

                           # dataframe columns
                           weighting='logFC', 
                           cluster='group',
                           gene='feature', 

                           # other options
                           show_num=30, 
                           value='avgExpr', 
                           AP=NULL, 
                           organism_db=NULL
                           ){
        AP <- return_aes_param (AP)
        if (is.null(markers)){
                markers <- find_DE_genes (x, save_dir, group.by=group.by, label=label)
        }
        if (show_num != 'all'){
        markers %>% dplyr::group_by (!!as.symbol (cluster) ) %>%
                dplyr::top_n (n=show_num, wt=!!as.symbol (weighting) ) %>%
                dplyr::ungroup () %>% as.data.frame () -> top_markers
        }else{top_markers <- markers}

        cell_types <- unique (top_markers [, cluster])
        cell_types <- sort (partial_relevel (cell_types, AP$cell_order))
        save_name <- paste ('DE_', group.by, '_', label, '.xlsx', sep='')
        save_path <- paste (save_dir, save_name, sep='/')

        if (save_format == 'excel'){
                save_DE_excel (top_markers, cell_types, weighting, cluster, save_path)
        }else{
                save_path <- gsub ('.xlsx', '_sim.csv', save_path)
                save_DE_gene_table (top_markers, value, cluster, gene, save_path, AP, organism_db)
        }
}

#' @importFrom magrittr %>%
save_DE_excel <- function (top_markers, cell_types, weighting, cluster, save_path){
        for (i in 1:length (cell_types)){
                if (i==1) {append_sheet <- F}else{append_sheet <- T}
                top_markers [top_markers [, cluster ] == as.character (cell_types[i]),] %>%
                        dplyr::arrange (desc (!!as.symbol (weighting)  )) %>%
                        xlsx::write.xlsx (file=save_path, 
                        sheetName = as.character (cell_types [i]), append=append_sheet)
        }
}

#' @importFrom magrittr %>%
save_DE_gene_table <- function (top_markers, value, cluster, gene, save_path, AP, organism_db){
        data (TF, package='TBdev')
        top_markers %>% dplyr::select (dplyr::all_of (c(value, cluster, gene)) ) %>%
                magrittr::set_colnames(c('avgExpr', 'cluster', 'gene')) %>%
                dplyr::mutate (cluster =partial_relevel (cluster, AP$cell_order)) %>%
                reshape2::dcast (gene ~ cluster, mean, value.var='avgExpr') %>%
                dplyr::mutate (transcription_factor = gene %in% TF) -> final_df
        if (!is.null (organism_db)) {
                final_df$ensemble_ID <- AnnotationDbi::mapIds (organism_db, 
                        keys=as.character (final_df$gene), keytype='SYMBOL',
                        column='ENSEMBL')
        }
        utils::write.csv (final_df, save_path)
}


#' Save pairwise DE gene in excel files
#' 
#' @description As many xlsx files will be created as there are comparison
#' groups. Each file contains the comparison groups
#' @param ... arguments passing to `save_DE_genes`
#' @importFrom magrittr %>%
save_DE_pairwise <- function (DE_df, group.by, directory, ...){
        all_types <- DE_df$compare_group %>% unique () %>% as.character () 
        if (!dir.exists (directory)){dir.create (directory)}
        for (i in all_types){
                DE_df [DE_df$compare_group == i,] %>% 
                        dplyr::select (!compare_group) -> sel_DE
                save_DE_genes (x=NULL, save_dir=directory, group.by=group.by, 
                               markers=sel_DE, label=i,...)
        }
}

# --------
# Plotting
# --------

#'  Visualise DE genes
#'
#' @param x Seurat object
#' @param markers DE gene marker dataframe computed by
#' `FindAllMarkers`. Alternatively, one can parse a vector of marker
#' genes.
#' @param top_gene how many genes to show per group
#' @param by_group which group of features in which the top DE genes
#' are selected
#' @param show_group which group of features in which the expression of
#' top DE genes are plotted against
#' @param show_heatmap whether to show heatmap or violin plot
#' @export
plot_DE_genes <- function (x, markers, top_gene=1, by_group='cluster',
                           show_group='Type', show_heatmap=TRUE){
        if (class (markers) == "data.frame"){
                top_markers <- unique_DE_genes (markers, top_gene)
                feature_plot <- as.character (top_markers$gene)
        }else{feature_plot <- as.character (markers)}

        print (feature_plot)
        if (!show_heatmap){
                plot_data <- Seurat::VlnPlot(x, features = feature_plot , group.by=show_group, ncol
                        = length(feature_plot)%/%3)
        }else{
                plot_data <- Seurat::DoHeatmap(x, features = feature_plot, group.by=show_group) + 
                        Seurat::NoLegend()
        }
        return (plot_data)
}

#' Similar to Seurat VlnPlot, with more graphic controls
#'
#' @param x a Seurat object
#' @param features genes to plot, which will be split over by facets
#' @param group.by by which cell feature
#' @param num_col in how many columns the facets are arranged
#' @param ... additional arguments for `custom_tick`
#' @importFrom magrittr %>%
#' @importFrom ggplot2 aes aes_string
#' @export
seurat_violin <- function (x, features, group.by, assay='RNA',
                           slot_data='data', num_col=NULL, free_xy='fixed',
                           AP=NULL, box_plot=T, ...){
        AP <- return_aes_param (AP)
        x %>% Seurat::FetchData (vars=c(features, group.by)) %>%
                magrittr::set_colnames (c(features, 'feature')) %>%
                tidyr::gather ( 'gene', 'expr_val', -feature ) %>%
                tidyr::drop_na () %>%
                dplyr::mutate ( gene = factor(gene, levels=features) ) -> plot_data

        ggplot2::ggplot (plot_data, aes ( x= feature, y=expr_val, fill=feature) ) +
                ggplot2::facet_wrap (.~gene, ncol=num_col, scales=free_xy) +
                ggplot2::labs (fill='')+ ggplot2::xlab ('')+
                ggplot2::ylab (expression (paste (italic ('mRNA levels')))) +
                theme_TB ('dotplot', feature_vec=x@meta.data[, group.by],
                          color_fill=T, aes_param=AP, rotation=90)+
                ggplot2::guides( fill= ggplot2::guide_legend(override.aes = list(alpha=1)))+
                custom_tick (min_prec=1,...) -> plot_ob
        if (box_plot){plot_ob <- plot_ob + ggplot2::geom_boxplot()
        }else{plot_ob <- plot_ob +ggplot2::geom_jitter (shape=AP$normal_shape, 
                                size=AP$pointsize, height=0, color=AP$point_edge_color, 
                                stroke=AP$edge_stroke)}
        return (plot_ob)
}

#' violin plot on genes that contribute to latent dimensions
strong_gene_violin <- function (x, group.by, dim_num=1, DR='pca', assay='RNA',
                                slot_data='data'){
        feature_load <- x@reductions[[DR]]@feature.loadings
        data (TF, package='TBdev')
        strong_gene <- data.frame (feature_load [rownames (feature_load) %in% TF [,1], dim_num])
        abs (strong_gene) %>% dplyr::top_n (20) -> gene_plot
        vln_plot <- seurat_violin (x, rownames (gene_plot), 'species_type', assay, slot_data)
        return (vln_plot + ggplot2::ggtitle (colnames (feature_load)[dim_num] ) )
}

#' Deseq2 style Volcano plot
#'
#' @param markers a dataframe, result generated from `find_DE_genes`
#' @param group1 comparison group1
#' @param group2 comparison group2
#' @param group1_col a column in `markers` where `group1` information is stored
#' @param group2_col a column in `markers` where `group2` information is stored
#' @param label_genes which genes to show in the volcano plot
#' @param weighting the x axis, i.e. log fold change
#' @param pval the p value for the y axis
#' @param gene where the gene names are stored
#' @param weight_thres above which value in `weighting` are the points colored
#' according to their comparison groups
#' @param AP aesthetic parameter for controlling the plot
#' @param ... arguments for `enrich_arrow`, including `nudge_x`, `nudge_y` and
#' `band_thickness`, `length_ratio`. For description of their roles, see
#' `enrich_bar`
#' @return a ggplot object
#' @importFrom magrittr %>%
#' @importFrom ggplot2 aes aes_string
#' @export
seurat_volcano <- function (markers, group1, group2, group1_col = 'group',
                            group2_col = 'compare_group', label_genes=NULL,
                            weighting='logFC', pval='padj', gene='feature',
                            weight_thres=0.25, show_gene_num=30, logy=T, AP=NULL, ...){
        AP <- return_aes_param (AP)
        # select appropriate data 
        markers %>% dplyr::filter (!!as.symbol (group1_col) == group1 ) %>%
                dplyr::filter (!!as.symbol (group2_col) == group2 ) -> sel_mark
        if (logy){
                sel_mark %>% dplyr::mutate (logp = - log10 (!!as.symbol (pval) )) -> plot_data
                y_lab <- expression ('-log[10] P')
        }else{
                sel_mark %>% dplyr::mutate (logp = !!as.symbol (pval)) -> plot_data
                y_lab <- pval
        }

        # color the significant points
        plot_data$groupby <- 'unknown'
        plot_data$groupby [plot_data [, weighting] > weight_thres] <- group1
        plot_data$groupby [plot_data [, weighting] < -weight_thres] <- group2

        # remove the Inf logp values
        plot_data %>% dplyr::filter (logp!=Inf) %>% dplyr::pull (logp) %>% max () -> max_p
        plot_data$logp [plot_data$logp==Inf] <- max_p

        # basic volcano plot
        ggplot2::ggplot (plot_data, aes_string (x=weighting, y='logp') )+
                ggplot2::geom_point (aes (fill=groupby), shape=AP$normal_shape, 
                                     color=AP$point_edge_color, show.legend=F,
                                     size=AP$pointsize) +
                # construct the threshold line
                ggplot2::geom_vline (xintercept=c(-weight_thres, weight_thres),
                                     linetype='dashed') +
                # add arrow sign
                enrich_arrow (plot_data, AP, yval='logp', xval=weighting,
                              group1='group', char_y=F, ...)+
                ggplot2::ylab (y_lab) -> plot_ob

        # add label genes
        if (is.null (label_genes) & !is.null (show_gene_num) ){
                sel_mark %>% dplyr::mutate (pos_neg = ifelse (!!as.symbol (weighting) > 0, 'pos', 'neg') ) %>%
                dplyr::group_by (pos_neg) %>% 
                dplyr::slice_max (abs (!!as.symbol (weighting)), n=show_gene_num/2) %>% 
                dplyr::pull (!!as.symbol(gene)) -> label_genes
        }
        if (!is.null (label_genes)){
                plot_data %>% dplyr::filter (!!as.symbol (gene) %in% label_genes) -> text_df
                if (!is.null (names (label_genes) )){
                        text_df$color_by <- names (label_genes) [match (text_df[, gene], label_genes)]
                        text_df$color_by <- partial_relevel (text_df$color_by, AP$cell_order)
                }else{
                        data.table::fifelse (text_df [, weighting] > 0, group1,
                                             group2) -> text_df$color_by 
                }
                aes_arg_list <- list (label=gene, color='color_by')
                aes_arg <- do.call (aes_string, aes_arg_list)
                plot_ob <- plot_ob + 
                        ggrepel::geom_text_repel(aes_arg, data=text_df, 
                                size=AP$point_fontsize, family=AP$font_fam, show.legend=F)
                feature_vec <- text_df$color
        }else{feature_vec <- NULL}
        plot_ob + theme_TB ('dotplot', rotation=0, aes_param=AP, color_fill=F,
                            feature_vec=feature_vec)+
        add_custom_color (feature_vec=plot_data$groupby, aes_param=AP, color_fill=T)
}

#' Gene expression in two categories
#'
#' @description Calculate the mean expression of each gene in two categories
#' @param x a Seurat object
#' @param compare.by which features to be compared
#' @param group.by perform the analysis separately for different groups of
#' cells
#' @return a dataframe
#' @importFrom magrittr %>%
gene_gene_mat <- function (x, compare.by, group.by, assay='RNA', slot_data='data'){
        datExpr <- as.matrix (Seurat::GetAssayData (x, assay=assay, slot=slot_data) )
        group_compare <- paste ( x@meta.data [, compare.by], x@meta.data [, group.by], sep='__')
        unique_group <- unique (group_compare)
        mean_expr <- lapply ( as.list (unique_group), function (inp) { 
                                     rowMeans (datExpr [, group_compare == inp, drop=F] ) } )
        mean_expr <- do.call (cbind, mean_expr)
        colnames (mean_expr) <- unique_group

        t (mean_expr) %>% data.frame () %>%
                tibble::add_column (compare.group=unique_group) %>%
                tidyr::separate (compare.group, c('compare_item', 'group_item'), sep='__') %>%
                tidyr::gather ('gene', 'expr_val', -compare_item, -group_item) %>%
                tidyr::spread (compare_item, expr_val) %>%
                tidyr::drop_na () 
}
