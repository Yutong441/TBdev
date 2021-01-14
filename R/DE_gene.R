# devtools::install_github('immunogenomics/presto')

# ----------------
# DE gene analysis
# ----------------

#' Find DE genes across a particular feature
#' 
#' @param x a Seurat object
#' @param feature a column in the meta data
#' @param directory where to save the computed markers
#' @param load_data whether to load existing marker data
#' @export
find_DE_genes <- function (x, directory, feature='seurat_clusters', filename=
                           NULL, load_data=TRUE, label='all'){

        if (is.null (filename)){
                filename <- paste ('marker_cluster_', feature, '_', label, '.csv', sep='')
        }
        file_name <- paste (directory, filename, sep='/')

        if (!load_data | !file.exists (file_name) ){
                Idents (x) <- x@meta.data [, feature]
                # Don't use the `FindMarkers` function in Seurat which is too
                # slow
                markers <- presto::wilcoxauc (x, feature, assay='data'  )
                utils::write.csv (markers, file_name)
        }else{
                markers <- utils::read.csv (file_name, row.names=1)
        }
        return (markers)
}

find_DE_genes_selected <- function (x, save_dir, feature, compare_against, label){
        x$compare <- as.character (x@meta.data [, feature])
        x$compare [compare_against] <- 'rest'
        Idents (x) <- x$compare

        marker_list <- list ()
        for (ident_1 in unique (x$compare) ){
                if (ident_1 != 'rest'){
                        print (paste ('processing', ident_1) )
                        marker_list[[ident_1]] <- Seurat::FindMarkers (x, ident.1=ident_1,
                                                 ident.2='rest', assay='RNA',
                                                 slot='data')
                        marker_list[[ident_1]]$cluster <- ident_1
                }
        }
        names (marker_list) <- rep('id', length (marker_list) )
        marker_all <- do.call (rbind, marker_list)
        rownames (marker_all) <- gsub ('id.', '', rownames (marker_all))
        write.csv (marker_all, paste (save_dir, paste ('marker_cluster_', 
                                feature, '', label, '.csv', sep=''), sep='/') )
} 

#' Select the top DE genes from each cluster
#' If 1 DE gene appears in 2 clusters, the cluster with the highest avg_logFC
#' will be chosen
#' NB: only upregulated genes will be selected
#'
#' @param DE_genes a dataframe
#' @param top_number how many DE genes per cluster/group to choose
#' @param gene which column has the gene symbols
#' @param weighting which column has the gene rankings
#' @param cluster which column has the cluster names/numbers
#' @importFrom magrittr %>%
#' @export
unique_DE_genes <- function (DE_genes, top_number, gene='feature',
                             weighting='logFC', cluster='group'){
        DE_genes %>%
                dplyr::group_by (!!as.symbol (cluster) ) %>%
                dplyr::top_n (n=top_number, wt=!!as.symbol (weighting) ) %>%
                dplyr::ungroup () %>%
                dplyr::arrange (desc (!!as.symbol (weighting) ) ) %>%
                dplyr::distinct (!!as.symbol (gene), .keep_all=T ) -> top_markers
        return (top_markers)
}

#' Save Differentially Expressed Genes
#'
#' @param x a Seurat object
#' @importFrom magrittr %>%
#' @export
save_DE_genes <- function (x, save_dir, feature='revised', label='all', 
                           show_num=30, weight='logFC', cluster='group'){
        markers <- find_DE_genes (x, save_dir, feature=feature, label=label)
        markers %>% group_by (!!as.symbol (cluster) ) %>%
                top_n (n=show_num, wt=!!as.symbol (weight) ) %>%
                ungroup () %>% as.data.frame () -> top_markers

        cell_types <- unique (top_markers [, cluster])
        save_name <- paste ('Sum_table_DE_', feature, '_', label, '.xlsx', sep='')
        for (i in 1:length (cell_types)){
                if (i==1) {append_sheet <- F}else{append_sheet <- T}
                top_markers [top_markers [, cluster ] == as.character (cell_types[i]),] %>%
                        dplyr::arrange (desc (!!as.symbol (weight)  )) %>%
                        xlsx::write.xlsx (file=paste (save_dir, save_name, sep='/'), 
                        sheetName = as.character (cell_types [i]), append=append_sheet)
        }
}

list_DE_genes <- function (x, directory, feature, top_number=8, label='all',
                           gene='feature', cluster='group'){
        DE_genes <- find_DE_genes (x, directory, feature, label=label)
        top_markers <- as.data.frame (unique_DE_genes (DE_genes, top_number))
        DE_genes <- as.vector (top_markers [, gene])
        names (DE_genes) <- top_markers[,cluster]
        return (DE_genes)
}

#' Merge DE genes and known lineage markers
#'
#' @export
DE_lineage_genes <- function (x, directory, top_number=3, label='all',
                              feature='revised', lineage_markers=NULL){
        if (is.null (lineage_markers)){data (lineage_markers) }
        DE_genes <- list_DE_genes (x, directory, feature=feature,
                                   top_number=top_number, label=label)
        all_genes <- c(lineage_markers,  DE_genes)
        all_genes <- all_genes [ gtools::mixedorder (names (all_genes)) ]
        return (all_genes)
}

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
#' @importFrom magrittr %>%
#' @importFrom ggplot2 aes aes_string
#' @export
seurat_violin <- function (x, features, group.by, assay='RNA',
                           slot_data='data', num_col=3, AP=NULL){
        AP <- return_aes_param (AP)
        datExpr <- as.matrix ( Seurat::GetAssayData (x, assay=assay, slot=slot_data)   ) 
        datExpr <- datExpr [rownames (datExpr) %in% features, ]
        print ('violin plot')
        plot_ob <- t (datExpr) %>% as.data.frame () %>%
                tibble::add_column (feature = x@meta.data [, group.by]) %>%
                tidyr::gather ( 'gene', 'expr_val', -feature ) %>%
                dplyr::mutate ( gene = factor(gene, levels=features) ) -> plot_data

        plot_ob <- ggplot2::ggplot (plot_data, aes ( x= feature, y=expr_val, fill=feature) ) +
                #geom_violin (show.legend=F) + 
                ggplot2::geom_jitter (shape=21, color='white', size=AP$pointsize) +
                ggplot2::facet_wrap (.~gene, ncol=num_col) +
                ggplot2::theme (axis.text.x = element_text (angle=90) ) +
                ggplot2::labs (fill='')+ ggplot2::xlab ('')+
                theme_TB ('dotplot', feature_vec=x@meta.data[, group.by], color_fill=T, AP=AP)+
                custom_tick (plot_data$expr_val) 
        return (plot_ob)
}

#' violin plot on genes that contribute to latent dimensions
strong_gene_violin <- function (x, group.by, dim_num=1, DR='pca', assay='RNA', data_slot='data'){
        feature_load <- x@reductions[[DR]]@feature.loadings
        data (TF)
        strong_gene <- data.frame (feature_load [rownames (feature_load) %in% TF [,1], dim_num])
        abs (strong_gene) %>% top_n (20) -> gene_plot
        vln_plot <- seurat_violin (x, rownames (gene_plot), 'species_type', assay, data_slot)
        return (vln_plot + ggplot2::ggtitle (colnames (feature_load)[dim_num] ) )
}

#' Plot the gene 
#'
#' @param x a Seurat object
#' @param compare.by which features to be compared
#' @param compare_items which 2 features in `compare.by` to be compared
#' @param group.by perform the analysis separately for different groups of
#' cells
#' @importFrom magrittr %>%
#' @export
gene_gene_mat <- function (x, compare.by, compare_items, group.by, assay='RNA', data_slot='data'){
        datExpr <- as.matrix (GetAssayData (x, assay=assay, slot=data_slot) )
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
