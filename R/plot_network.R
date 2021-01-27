#' Plot gene-gene correlation network in a WGCNA cluster
#' 
#' @param x a Seurat object
#' @param color_row a vector of gene names
#' @param cell_type the cell type in which correlation is computed. By default,
#' all cell types will be used
#' @param feature which metadata entry the cell type information is stored.
#' Ignore it if `cell_type='all'`
#' @param scale_edge_weight adjust the weight of the lines
#' @param scale_node_size adjust the relative node size
#' @param threshold only gene-gene pair with correlation beyond a threshold
#' will be selected
#' @param return_igraph whether igraph object will be returned or whether it
#' will be plotted
#' @importFrom igraph E "E<-"
#' @importFrom magrittr %>%
#' @export
plot_WGCNA_net <- function (x, color_row, celltype='all', feature='revised',
                            scale_edge_weight=0.2, scale_node_size=0.2,
                            threshold=0.5, return_igraph=F){
        if (return_igraph){scale_edge_weight <- 1}
        datExpr <- subset_celltype (x[color_row,], 'all', celltype, featurue)
        sim_mat <- WGCNA::cor (datExpr, method="pearson")
        na_field <- apply (sim_mat, 1, function(x){mean (is.na (x)) }) == 1
        sim_mat <- sim_mat [!na_field, !na_field]

        sim_mat %>% data.frame () %>%
                tibble::add_column (gene1 = rownames (sim_mat) ) %>%
                tidyr::gather ('gene2', 'val', -gene1) %>%
                dplyr::filter (val > threshold) -> select_sim
        select_sim <- select_sim [select_sim$gene1 != select_sim$gene2,]
        sel_genes <- unique (select_sim$gene1)
        sim_mat <- sim_mat [sel_genes, sel_genes]

        network <- igraph::graph_from_adjacency_matrix(sim_mat , 
                        weighted=T, mode="undirected", diag=F)
        E (network)$width <- abs (E (network)$weight/scale_edge_weight)
        E (network)$weight_sign <- sign (E (network)$weight)
        E (network)$weight <- abs (E(network)$weight)

        avg_exp <- colMeans (datExpr)
        avg_exp <- avg_exp [ match (rownames (sim_mat), names (avg_exp) ) ]
        if (!return_igraph){
                igraph::plot.igraph (network, vertex.size=as.matrix (
                                                avg_exp)/scale_node_size)
        }else{
                V (network)$size <- as.matrix (avg_exp)
                return (network)
        }
}

#' Plot the gene-gene interaction network for all WGCNA clusters
#'
#' @description This is basically repeating the `plot_WGCNA_net` function and
#' save every plot in a single pdf file. 
#' @param x a Seurat object
#' @param color_row a list of vector of gene names
#' @param save_dir where to save the plot
#' @export
plot_all_WGCNA_nets <- function (x, all_groups, save_dir, ...){
        grDevices::pdf (save_dir)
        lapply (all_groups, function (one_group){
                plot_WGCNA_net (x, one_group, ...) 
                        })
        grDevices::dev.off()
}

#' Customise igraph objects using ggraph
#'
#' @param igraph_net an igraph object
#' @param manual_xy the x and y coordinates for the plots
#' @param order_leftright order the points from left to right according to the
#' a feature of the vertices, e.g. pseudotime. If this is set to be NULL, no
#' reordering will occur unless `manual_xy` is not NULL. If `order_leftright`
#' is not NULL, the setting will override that in `manual_xy`.
#' @export
custom_net <- function (igraph_net, AP=NULL, limits=NA, ranges=c(0,15),
                        manual_xy=NULL, order_leftright='pseudotime',
                        ascend_order=T){
        AP <- return_aes_param (AP)
        if (is.na (limits[1])){limits <- range (V(igraph_net)$size) }
        if (!is.null(order_leftright)){
                xy <- graphlayouts::layout_with_stress(igraph_net)
                xy_new <- xy [order (xy[,1]),]
                old_genes <- attr (V(graph_net), 'names')
                new_genes <- old_genes [order (V(graph_net)$pseudotime, 
                                               decreasing=!ascend_order)]
                manual_xy <- xy_new [match (old_genes, new_genes),]
        }
        if (is.null(manual_xy)){
                gnet <- ggraph::ggraph (igraph_net, layout='kk') 
        }else{
                gnet <- ggraph::ggraph (igraph_net, 'manual', x=manual_xy[,1], y=manual_xy[,2]) 
        }
        gnet+   ggraph::geom_edge_link (ggplot2::aes(edge_width=width), edge_color='gray') +
                ggraph::geom_node_point (ggplot2::aes (size=size, fill=pseudotime), shape=21)+
                ggraph::geom_node_text (ggplot2::aes(label=name), vjust='bottom', 
                                        repel=T, size=AP$point_fontsize) +
                TBdev::theme_TB ('no_arrow')+ ggplot2::scale_fill_viridis_c ()+
                ggplot2::scale_size (range=ranges, limits=limits)+
                ggplot2::labs (size='norm count')
}

#' @export
#' @importFrom igraph E "E<-"
#' @importFrom magrittr %>%
add_pseudotime_to_net <- function (x, igraph_net, AP=NULL, assay='RNA',
                                 slot_data='data',
                                 time_ref='assigned_cluster', PT_col='MGP_PT'){
        sel_genes <-attr (V(igraph_net), 'names')
        Seurat::GetAssayData (x[sel_genes, ], assay=assay, slot=slot_data) %>%
                as.matrix () %>% t() %>% data.frame () %>%
                tibble::add_column (Type=x@meta.data[, time_ref]) -> pt_dat

        pt_dat %>% tidyr::gather ('gene', 'expre', -Type) %>%
                dplyr::group_by (gene, Type) %>% 
                dplyr::summarise (mean_expr=mean(expre)) %>%
                dplyr::ungroup () %>% dplyr::group_by (gene) %>% 
                dplyr::slice_max(mean_expr, n=1) %>% 
                data.frame () -> gene_type_df

        x@meta.data %>% dplyr::select (dplyr::all_of (c(PT_col, time_ref))) %>% 
                tidyr::drop_na () %>%
                magrittr::set_colnames (c('PT', 'Type')) %>%
                dplyr::group_by (Type) %>% dplyr::summarise (mean_pt = mean (PT)) %>% 
                data.frame () -> type_time_df
        gene_type_df$PT <- type_time_df$mean_pt [match (gene_type_df$Type, 
                                                        type_time_df$Type)]
        V(igraph_net)$pseudotime <- gene_type_df$PT [match (attr(V(igraph_net), 
                                                'names'), gene_type_df$gene)]
        return (igraph_net)
}

#' Add the mean expression level of a cell type to igraph
#'
#' @param markers a dataframe generated by `find_DE_genes`
#' @importFrom igraph E "E<-"
#' @importFrom magrittr %>%
#' @export
change_expre_level <- function (markers, igraph_net, celltype, group.by='group',
                                assay='RNA', slot_data='data',
                                expr_col='avgExpr', gene_col='feature'){
        sel_genes <-attr (V(igraph_net), 'names')
        markers %>% dplyr::filter (!!as.symbol(group.by) == celltype) %>%
                dplyr::filter (!!as.symbol(gene_col) %in% sel_genes) %>%
                dplyr::select (dplyr::all_of (c(gene_col, expr_col))) %>% 
                tibble::deframe () -> gene_expr
        V(igraph_net)$size <- gene_expr [match (sel_genes, names (gene_expr) )]
        return (igraph_net)
}

#' @export
custom_net_cells <- function (markers, igraph_net, celltypes, limits=NA,
                              AP=NULL,...){
        graph_list <- list()
        for (i in 1:length(celltypes)){
                one_graph <- change_expre_level (markers, igraph_net, celltypes[i])
                graph_list[[i]] <- custom_net (one_graph, limits=limits, ...)
        }
        return (ggpubr::ggarrange (plotlist=graph_list, labels=celltypes) )
}
