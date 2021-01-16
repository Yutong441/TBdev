# carry out over-representation test

non_default_theme <- function (plot_ob, default_theme, AP, color_fill=F){
        if (!default_theme){
                col_vec <- AP$heatmap_color
                plot_ob <- plot_ob +
                        ggplot2::scale_color_continuous (high=col_vec[1], low=col_vec[3]) +
                        theme_dotplot (rotation=0, aes_param=AP)
                if (color_fill){plot_ob <- plot_ob + 
                        ggplot2::scale_fill_continuous (high=col_vec[1], low=col_vec[3]) }
                return (plot_ob)
        }
}

#' @importFrom ggplot2 aes
#' @noRd
gg_dot <- function (xx, showCategory){
        AP <- return_aes_param (NULL)
        plot_data <- xx@compareClusterResult
        plot_data$Cluster <- partial_relevel (plot_data$Cluster, AP$cell_order)

        plot_data %>% 
                dplyr::arrange (p.adjust) %>%
                dplyr::slice_min (p.adjust, n=showCategory) %>%
                tidyr::separate (GeneRatio, c('numer', 'denom'), sep='/') %>%
                dplyr::mutate_at (c('numer', 'denom'), as.numeric ) %>%
                dplyr::mutate (GeneRatio = numer/denom) %>%
                dplyr::group_by (Cluster) %>%
                dplyr::arrange (dplyr::desc (GeneRatio)) %>%
                dplyr::mutate (Description = factor (Description, levels = Description ) ) %>%
                ggplot2::ggplot ( aes (x=Cluster, y=stats::reorder (Description, dplyr::desc (Description)), 
                              color=p.adjust, size=GeneRatio) ) +
                ggplot2::geom_point ()+ ggplot2::ylab('Description')
}


#' Over-representation test of genes in each cluster
#'
#' @param markers a dataframe of DE genes from Seurat; alternatively, a list of
#' DE genes may be used
#' @param GO_data result from go_data
#' @param organism_db which organism database, for example, org.Hs.eg.db
#' @param organism_name which organism for KEGG and Reactome
#' @param logFC_term which column in `markers` contains the information for
#' ranking gene expression
#' @param logFC_thres threshold of `logFC_term` beyond which a gene is selected
#' for over-representation test
#' @param enrich_area either KEGG or GO or reactome
#' @export
compare_cluster_enrichment <- function (markers, GO_data, organism_db, organism_name='human', 
                                        logFC_term = 'logFC',
                                        log_FC_thres=0.25, enrich_area='KEGG'){
        if (class (markers) == 'data.frame'){
                markers %>% filter ( !!as.symbol (logFC_term) > log_FC_thres ) -> markers
                cell_type <- unique (markers$group)
                markers$entrez <- AnnotationDbi::mapIds(organism_db, as.character (
                                                markers$feature), 'ENTREZID', 'SYMBOL')
                gene_list <- lapply ( as.list (cell_type), function (x){
                                             markers$entrez [markers$group == x]} )
                names (gene_list) <- cell_type
        }else{
                gene_list <- lapply (markers, function (x) {AnnotationDbi::mapIds (
                                                organism_db, as.character (x), 
                                                'ENTREZID', 'SYMBOL') %>% unlist ()}) 
                names (gene_list) <- names (markers)
        }
        if (enrich_area == 'GO'){
                kk <- clusterProfiler::compareCluster(gene_list, 
                                     fun=paste ('enrich', enrich_area, sep=''), 
                                     OrgDb=organism_db, pvalueCutoff=0.05)
        }else if (enrich_area == 'reactome'){
                kk <- compare_reactome (gene_list)
        }else{
                kegg_name <- get_kegg (organism_name)
                kk <- clusterProfiler::compareCluster(gene_list, 
                                     fun=paste ('enrich', enrich_area, sep=''), 
                                     organism=kegg_name, pvalueCutoff=0.05)
        }
        return (enrichplot::pairwise_termsim(kk, method="JC", semData = GO_data))
}

#' Compare reactome across multiple cell types
#'
#' The `compareCluster` function does not work for 'enrichPathway' in my
#' computer. I need to write a similar function.
compare_reactome <- function (gene_list, cutoff=0.05){
        ra <- lapply (gene_list, function (x) {ReactomePA::enrichPathway (x, organism='human', 
                                                                          pvalueCutoff=cutoff )}) 
        ra_df <- lapply (names (gene_list), function (x) {ra [[x]]@result$Cluster <- x; 
                         return (ra[[x]]@result) } )
        ra_df <- do.call (rbind, ra_df)
        # the pvalueCutoff argument in the original function does not seem to work
        ra_df %>% dplyr::filter (p.adjust < cutoff) -> ra_df
        ra_obj <- methods::new ('compareClusterResult', compareClusterResult=ra_df)
        return (ra_obj)
}

#' Visualise results from over-representation test
#'
#' @param xx result from `compareCluster` or just a dataframe
#' @param show_graph whether to show the terms in piechart format ('emap') or
#' dotplot ('dotplot') or ridgeplot ('ridgeplot')
#' @param default_theme whether to use ggplot default theme
#' @param feature_vec to manually assign color
#' @param enrich_area either 'GO', 'KEGG' or 'reactome'. Only important for
#' saving the figure
#' @param clean_results whether to remove unrelated and disease related terms
#' @param simplification whether to simplify the GO/KEGG/reactome terms
#' @param sim_dict a data frame with 2 columns: 'ori' for the original terms,
#' 'sub' for the strings that will replace the original terms. The default is a
#' a limited data frame built into this package.
#' @param append_default_dict append default dictionary to simplify the terms
#' @description need to set seed to ensure reproducibility
#' @export
display_cluster_enrichment <- function (xx, show_graph='emap',
                                        feature_vec=NULL, default_theme=F,
                                        save_dir=NULL, enrich_area='GO',
                                        show_num=NULL, clean_results=T,
                                        simplification=T, AP=NULL,
                                        sim_dict=NULL,
                                        append_default_dict=T,...){
        AP <- return_aes_param (AP)
        sim_dict <- append_default_dictionary (sim_dict, append_default_dict)
        if (clean_results){ 
                ori_terms <- nrow (xx@compareClusterResult)
                xx <- clean_terms (xx, AP) 
                new_terms <- nrow (xx@compareClusterResult)
                print ( paste ('removing', ori_terms-new_terms,'non development related terms' ))
        }
        if (simplification){ rename_vec <- simplify_terms (xx, sim_dict)
        }else{ rename_vec <- NULL }
        #if (!is.null (feature_vec)){default_theme<-F}
        if (show_graph == 'emap'){
                if (is.null (show_num)){show_num=150}
                cluster_plot <- sim_emap (xx, pie='count', rename_vec=rename_vec,
                                        showCategory=show_num, layout='kk', legend_n=1, 
                                        AP=AP,...) 
                cluster_plot <- cluster_plot + theme_dim_red (AP, color_fill=T)
                if (is.null(feature_vec)){feature_vec <- xx@compareClusterResult$Cluster}
                cluster_plot <- cluster_plot + add_custom_color (feature_vec, AP, color_fill=T)
        }else if (show_graph == 'dotplot') { 
                if (is.null (show_num)){show_num=40}
                cluster_plot <- gg_dot (xx, showCategory=show_num)
                cluster_plot <- non_default_theme (cluster_plot, default_theme, color_fill=T, AP=AP)
        }else if (show_graph == 'ridgeplot') {
                if (is.null (show_num)){show_num=30}
                cluster_plot <- ridgeplot (xx, showCategory=show_num)
                cluster_plot <- non_default_theme (cluster_plot, default_theme, color_fill=T, AP=AP)
        }else{
                cluster_plot <- gseaplot (xx, geneSetID=1)
                cluster_plot <- non_default_theme (cluster_plot, default_theme, AP=AP)
        }
        if (is.null (save_dir) ){return (cluster_plot)
        }else{
                cluster_plot
                ggplot2::ggsave ( paste (save_dir, paste ('GSEA_', enrich_area, '_', 
                        show_graph, '.pdf', sep=''), sep='/'), width=16, height=8)
        }
}
