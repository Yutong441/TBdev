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

#' @importFrom igraph V
hide_node <- function (igraph_net, manual_xy=NULL, hide_node_thres=NULL){
        if (is.null(manual_xy)){
                manual_xy <- graphlayouts::layout_with_stress(igraph_net)
        }
        if (!is.null(hide_node_thres)){
                hide_vec <- V(igraph_net)$size < hide_node_thres
                node_name <- attr (V(igraph_net), 'names')
                manual_xy_hide <- manual_xy [!hide_vec,]
                igraph_hide <- igraph::delete_vertices (igraph_net, node_name[hide_vec])
                size_vec <- V(igraph_net)$size
                names (size_vec) <- node_name
                return (list (igraph_hide, manual_xy_hide))
        }else{
                return (list (igraph_net, manual_xy))
        }
}

order_net_leftright <- function (igraph_net, ascend_order){
        xy <- graphlayouts::layout_with_stress(igraph_net)
        xy_new <- xy [order (xy[,1]),]
        old_genes <- attr (V(igraph_net), 'names')
        new_genes <- old_genes [order (V(igraph_net)$pseudotime, 
                                       decreasing=!ascend_order)]
        return ( xy_new [match (old_genes, new_genes),])
}

#' Customise igraph objects using ggraph
#'
#' @param igraph_net an igraph object
#' @param manual_xy the x and y coordinates for the plots
#' @param order_leftright order the points from left to right according to the
#' a feature of the vertices, e.g. pseudotime. If this is set to be NULL, no
#' reordering will occur unless `manual_xy` is not NULL. If `order_leftright`
#' is not NULL, the setting will override that in `manual_xy`.
#' @param ascend_order order the pseudotime in the ascending order, i.e. the
#' lowest value will be colored purple, highest yellow according to the viridis
#' scale.
#' @param coloring a vector with the same length as the number of nodes that
#' provide the colors
#' @param nudge_ratio by how much the labels are radiated from the center
#' @param hide_node_thres below which node size are the nodes and the
#' associated connections hidden.
#' @return a ggplot object
#' @export
custom_net <- function (igraph_net, AP=NULL, limits=NA, ranges=c(0,15),
                        manual_xy=NULL, order_leftright='pseudotime',
                        ascend_order=T, coloring=NULL, nudge_ratio=0., 
                        hide_node_thres=NULL, plot_title=NULL){
        AP <- return_aes_param (AP)
        if (is.null(coloring)){coloring <- V(igraph_net)$pseudotime}
        if (is.na (limits[1])){limits <- range (V(igraph_net)$size) }
        if (is.null(manual_xy)){
                if (!is.null(order_leftright)){
                        manual_xy <- order_net_leftright (igraph_net, ascend_order)
                }else{manual_xy <- graphlayouts::layout_with_stress(igraph_net)}
        }
        com_list <- hide_node (igraph_net, manual_xy, hide_node_thres)
        igraph_net <- com_list[[1]]
        manual_xy <- com_list[[2]]

        ggraph::ggraph (igraph_net, 'manual', x=manual_xy[,1], y=manual_xy[,2]) -> gnet
        gnet +  ggraph::geom_edge_link (ggplot2::aes(edge_width=width), edge_color='gray') +
                ggraph::geom_node_point (ggplot2::aes (size=size, fill=pseudotime), 
                                         shape=21)+
                ggraph::geom_node_text (ggplot2::aes(label=name), 
                                        repel=F, size=AP$point_fontsize,
                                        nudge_x = gnet$data$x*nudge_ratio,
                                        nudge_y = gnet$data$y*nudge_ratio
                                        ) +
                ggplot2::scale_size (range=ranges, limits=limits)+
                ggplot2::labs (size='norm count')+
                ggplot2::guides(size=F) -> g_ob

        if (!is.null(plot_title)){g_ob <- g_ob + ggplot2::ggtitle (plot_title) }
        g_ob + TBdev::theme_TB ('no_arrow', feature_vec=coloring,
                                 color_fill=T, more_prec=3, aes_param=AP)
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

#' Delete those nodes with a small expression levels
#'
#' @param graph_ob an igraph object
#' @param size_thres the threshold of expression below which nodes are removed
#' @description This function removes nodes from graph completely. Thus the
#' underlying graph structure will change during ploting. If you wish to
#' preserve the graph structure, use `hide_node`
#' @return an igraph object
filter_expression <- function (graph_ob, size_thres){
        node_names <- attr (V(graph_ob), 'names')
        small_nodes <- V(graph_ob)$size < size_thres
        igraph::delete_vertices(graph_ob, node_names[small_nodes])
}

#' @importFrom igraph V
#' @export
custom_net_cells <- function (markers, igraph_net, celltypes, limits=NULL,
                              AP=NULL, return_sep=F, normalize_cell=NULL, ...){
        if (!is.null (normalize_cell)){
                igraph_net <- change_expre_level (markers, igraph_net, normalize_cell)
                normalize_vec <- V(igraph_net)$size
        }
        graph_list <- list()
        for (i in 1:length(celltypes)){
                one_graph <- change_expre_level (markers, igraph_net, celltypes[i])
                if (!is.null (normalize_cell)) {
                        V(one_graph)$size <- V(one_graph)$size/normalize_vec
                }
                graph_list[[i]] <- one_graph
        }
        xy <- order_net_leftright (igraph_net, ascend_order=T)
        if (is.null(limits)){limits <- range (do.call (c, lapply (graph_list, function (gx){V(gx)$size} ) )) }
        graph_list <- lapply (graph_list, function (x){custom_net (x, limits=limits, manual_xy=xy,...)})
        if (!return_sep){
                return (ggpubr::ggarrange (plotlist=graph_list, labels=celltypes) )
        }else{
                return (graph_list)
        }
}

#' Create an igraph object from seurat
#'
#' @description a convenient wrapper for `plot_WGCNA_net`, then
#' `add_pseudotime_to_net` and `change_expre_level` and `filter_expression`
#' @param ... arguments to `add_pseudotime_to_net`
custom_net_from_seurat <- function(x, genes, markers, celltype, thres=0.6,
                                   size_thres=0., ...){
        graph_net <- plot_WGCNA_net (x, genes, return_igraph=T, thres=thres)
        graph_net <- add_pseudotime_to_net (x, graph_net, ...)
        graph_net <- change_expre_level (markers, graph_net, celltype)
        graph_net <- filter_expression (graph_net, size_thres)
        return (graph_net)
}

#' Multiple network graphs
#'
#' @param x a Seurat object
#' @param size_thres filter out nodes with expression levels smaller than a
#' certain threshold
#' @param ... argments to pass to `custom_net`
#' @importFrom igraph V
#' @export
custom_net_diff_nets <- function (x, gene_list, markers, size_thres=0., ...){
        all_genes <- do.call(c,gene_list)
        limits <- range (x[all_genes,][['RNA']]@data)

        # plot net
        graph_list <- lapply (as.list(1:length(gene_list)), function(i){
               custom_net_from_seurat (x, gene_list[[i]], markers, names
                                       (gene_list)[i], size_thres=size_thres)
                              })

        pt_limits <- do.call (c, lapply (graph_list, function(gx){V(gx)$pseudotime}) )
        lapply (as.list (1:length(graph_list)), function (i){
                        custom_net (graph_list [[i]], limits=limits, 
                                    coloring=pt_limits, 
                                    plot_title=names (gene_list)[i],...)})
}
