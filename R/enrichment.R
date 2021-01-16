#' Obtain gene list from DE gene results
#'
#' @description This function obtains a list of genes, the input format to GSEA
#' or over-representation test, from the top DE genes. For parameter setting,
#' see `run_GSEA`
get_gene_list <- function (markers, cluster_num, organism_db, gene_col='feature',
                           cluster_col='group', FC='logFC'){
        geneList <- markers [markers [, cluster_col] == cluster_num, FC]
        names (geneList) <- markers [markers [, cluster_col]== cluster_num, gene_col]
        names (geneList) <- AnnotationDbi::mapIds(organism_db, 
                                names(geneList), 'ENTREZID', 'SYMBOL')
        geneList <- geneList [!duplicated (names (geneList) ) ]
        return (geneList)
}

#' GSEA
#'
#' @param markers a dataframe with at least 3 columns: one for gene names
#' (`gene_col`), one cluster (`cluster_col`) and one for ranking (`FC`). If
#' they are not 'feature', 'group' and 'logFC' respectively, they need to be
#' supplied as keyword arguments in ...
#' @param cluster_num which cluster to compute enrichment
#' @param organism_db which organism database, for example, org.Hs.eg.db
#' @param organism_name which organism for KEGG and Reactome
#' @param enrich_area either 'GO', 'KEGG' or 'reactome'
#' @export
run_GSEA <- function (markers, cluster_num, organism_db, organism_name
                      ='human', enrich_area='GO', ...){
        geneList <- get_gene_list (markers, cluster_num, organism_db, ...)
        print (paste ('analysing', length (geneList), 'genes' ))
        if (enrich_area == 'GO'){
                ego <- clusterProfiler::gseGO(
                              geneList     =  sort (geneList, decreasing=T),
                              OrgDb        = organism_db,
                              ont          = "ALL",
                              nPerm        = 1000,
                              minGSSize    = 100,
                              maxGSSize    = 500,
                              pvalueCutoff = 0.05,
                              verbose      = FALSE)
                return (ego)
        }

        if (enrich_area == 'KEGG'){
                kegg_organism <- get_kegg (organism_name)
                kk <- clusterProfiler::gseKEGG(
                               geneList     = sort (geneList, decreasing=T),
                               organism     = kegg_organism,
                               minGSSize    = 120,
                               pvalueCutoff = 0.05,
                               verbose      = FALSE)
                return (kk)
        }
        if (enrich_area == 'reactome'){
                ra <- clusterProfiler::gsePathway(
                               geneList     = sort (geneList, decreasing=T),
                               organism     = organism_name,
                               minGSSize    = 120,
                               pvalueCutoff = 0.05,
                               verbose      = FALSE)
                return (ra)
        }
}

#' Combine the GSEA results from all clusters
#'
#' @param markers a dataframe containing the fold changes of genes with a
#' column specified as `cluster_col` that contains the cluster information
#' @param organism_db which organism database, for example, org.Hs.eg.db
#' @param show_num how many GO/KEGG terms to show. This function will select
#' those terms whose mean absolute fold changes are high
#' @param enrich_area either 'GO', 'KEGG' or 'reactome'
#' @param save_path save the GSEA results for later usage. GSEA is quite
#' computationally intensive, especially considering that it will be run for
#' all clusters
#' @importFrom magrittr %>%
#' @export
run_GSEA_all_types <- function (markers, organism_db, organism_name='human',
                                show_num=50, cluster_col='group',
                                enrich_area='KEGG', save_path=NULL, AP=NULL,
                                ...){
        AP <- return_aes_param (AP)
        if (is.null (save_path)  ){proceed <- T
        }else{
                if (!file.exists (save_path)){proceed <- T
                }else{proceed <- F}
        }
        if (proceed){
                cell_types <- unique (markers [, cluster_col])
                # run GSEA for each cluster
                gsea_list <- lapply (as.list (cell_types), function (x){
                                     run_GSEA (markers, x, organism_db,
                                               enrich_area=enrich_area,
                                               organism_name=organism_name,
                                               cluster_col=cluster_col, ...)})
                # obtain the dataframe
                gsea_df_list <- lapply (gsea_list, ridge_one_type)
                names (gsea_df_list) <- cell_types
                print (' combine the data')
                for (i in 1:length (gsea_df_list)){ 
                        if (nrow (gsea_df_list[[i]]) > 0){
                                gsea_df_list[[i]]$cluster <- names (gsea_df_list)[i] 
                        }
                }
                print ('combine before writing')
                gsea_df <- do.call(rbind, gsea_df_list)
                if (!is.null (save_path)) {utils::write.csv (gsea_df, save_path)}
        }else{ 
                print ('reading from precomputed GSEA. ')
                print ('If you do not like this feature, set save_path=NULL')
                gsea_df <- utils::read.csv (save_path, row.names=1) 
        }
        gsea_df %>%
                dplyr::mutate (abs_change = abs (value)) %>%
                dplyr::group_by (category) %>%
                dplyr::summarise (mean_change = mean (abs_change) ) %>%
                dplyr::slice_max (mean_change, n=show_num) -> selected_terms

        # filter out irrelevant terms
        relevant_terms <- remove_terms (gsea_df$category, AP)
        gsea_df %>%
                dplyr::filter (category %in% selected_terms$category) %>%
                dplyr::filter (category %in% relevant_terms) -> gsea_df
        return (gsea_df)
}


all_GSEA_one_type <- function (markers, cluster_num, save_dir, label=NULL, ...){
        if (is.null(label)){label <- cluster_num}
        save_dir <- paste (  save_dir, label, sep='/' )
        if (!dir.exists (save_dir) ){dir.create (save_dir) }

        ego <- run_GSEA (markers, cluster_num, organism_db, ...)
        display_cluster_enrichment (ego, show_graph='dotplot', save_dir=save_dir)
        display_cluster_enrichment (ego, show_graph='ridgeplot', save_dir=save_dir)
        display_cluster_enrichment (ego, show_graph='gseaplot', save_dir=save_dir)
        rm (ego)

        ekk <- run_GSEA (markers, cluster_num, organism_db, enrich_area='KEGG')
        display_cluster_enrichment (ekk, show_graph='dotplot', enrich_area='KEGG', save_dir=save_dir)
        display_cluster_enrichment (ekk, show_graph='ridgeplot', enrich_area='KEGG', save_dir=save_dir)
        display_cluster_enrichment (ekk, show_graph='gseaplot', enrich_area='KEGG', save_dir=save_dir)
        rm (ekk)
}

get_all_paths <- function (save_dir,  cluster_num, all_path=NULL, path_data_dir='GO/pathways'){
        if (is.null(all_path)){data (Kegg_ID, package='TBdev'); all_path <- Kegg_ID
        }
        final_dir <- paste (save_dir, cluster_num, sep='/')
        if (!dir.exists (final_dir) ){dir.create (final_dir) }
        for ( i in 1:nrow(all_path) ){
                path_ID <- all_path$kegg_id[i]
                old_name <-  paste (path_ID, '.pathview.png', sep='') 
                new_name <- paste ('Path_', all_path$pathway[i], '.png', sep='')
                new_name <-  paste (final_dir, new_name, sep='/') 
                file.rename (from=old_name, to=new_name)
        }
}

#' Ridgeplot
#'
#' @param x an enrichment object 
#' @references
#' \url{https://rdrr.io/github/GuangchuangYu/enrichplot/src/R/ridgeplot.R}
ridge_one_type <- function(x, fill='p.adjust', core_enrichment = T, orderBy =
                           "NES", decreasing = F) {
        if (!fill %in% colnames(x@result)) { stop("'fill' variable not available...") }
        if (orderBy !=  'NES' && !orderBy %in% colnames(x@result)) {
                message('wrong orderBy parameter; set to default `orderBy = "NES"`')
                orderBy <- "NES"
        }
        if (core_enrichment) { gs2id <- DOSE::geneInCategory(x)
        } else { gs2id <- x@geneSets[x$ID] 
        }
        gs2val <- lapply(gs2id, function(id) {
                res <- x@geneList[id]
                res <- res[!is.na(res)]
        })
        nn <- names(gs2val)
        i <- match(nn, x$ID)
        nn <- x$Description[i]

        j <- order(x@result[[orderBy]][i], decreasing = decreasing)
        len <- sapply(gs2val, length)
        gs2val.df <- data.frame(category = rep(nn, times=len),
                                color = rep(x[i, fill], times=len),
                                value = unlist(gs2val)
        )
        colnames(gs2val.df)[2] <- fill
        gs2val.df$category <- factor(gs2val.df$category, levels=nn[j])
        return (gs2val.df)
}

#' Ridge plot for multiple categories
#'
#' @param xx a dataframe generated from run_GSEA_all_types
#' @param x_col the column for the x axis of ridge plot
#' @param y_col the column for the y axis of ridge plot
#' @param sort_by by which column in `xx` the order of the enriched terms are
#' shown in ridge plot
#' @param default_theme whether to use ggplot default_theme
#' @param color_col color by which column in `xx`
#' @param not_show_prop proportion of the data with extreme values of fold
#' changes that are not shown, in order to truncate the x axis and avoid white
#' spaces in the final plot.
#' @param term_length remove terms beyond a certain length
#' @param show_num how many enriched terms to show
#' @param simplification whether to simplify the GO/KEGG/Reactome terms
#' @param sim_dict a data frame with 2 columns: 'ori' for the original terms,
#' 'sub' for the strings that will replace the original terms. The default is a
#' a limited data frame built into this package.
#' @param amplify_font increase the font size of the y axis labels by a ratio,
#' which represent the enrichment terms
#' @param AP aesthetic parameters
#' @param sim_dict a data frame with 2 columns: 'ori' for the original terms,
#' 'sub' for the strings that will replace the original terms. The default is a
#' a limited data frame built into this package.
#' @param append_default_dict append default dictionary to simplify the terms
#' @importFrom stats quantile
#' @importFrom ggplot2 aes_string
#' @importFrom magrittr %>%
#' @author Yutong Chen
#' @export
ridge_all_types <- function (xx, x_col='value', y_col='category', sort_by=
                             'p.adjust', default_theme=F, color_col='cluster',
                     not_show_prop=0.01, term_length=60, show_num=NULL,
                     simplification=T, amplify_font=1.2, AP=NULL,
                     sim_dict=NULL, append_default_dict=T){
        AP <- return_aes_param (AP)
        sim_dict <- append_default_dictionary (sim_dict, append_default_dict)

        xx %>% dplyr::mutate (char_length = nchar (as.character (!!as.symbol ('category')) ) ) %>%
                dplyr::filter (char_length < term_length) %>%
                dplyr::arrange (!!as.symbol (sort_by) ) -> plot_data
        if (simplification) {plot_data <- simplify_gsea (plot_data, sim_dict)} # from 'clean_terms.R'
        plot_data [, y_col] <- factor (plot_data [, y_col], levels=unique (plot_data [, y_col]))
        if (!is.null (show_num)){
                show_terms <- levels (plot_data [, y_col]) [1:show_num]
                plot_data %>% dplyr::filter ( !!as.symbol (y_col) %in% show_terms ) -> plot_data
        }

        # so that the most significant terms are in the top of the ridge plot
        plot_data [, y_col] <- factor (plot_data [, y_col], levels=rev (unique (plot_data [, y_col])) )
        xmin <- quantile (xx [, x_col], not_show_prop)[1]
        xmax <- quantile (xx [, x_col], 1-not_show_prop)[1]

        ggplot2::ggplot(plot_data, aes_string(x=x_col, y=y_col, fill=color_col)) +
                ggridges::geom_density_ridges(alpha=AP$ridge_alpha) +
                ggplot2::xlab(NULL) + ggplot2::ylab(NULL) +
                ggplot2::scale_x_continuous (breaks= round (seq(xmin, xmax, length.out=3 )), 
                                             limits=c(xmin, xmax))  -> plot_xx
        if (!default_theme){
                plot_xx <- plot_xx + theme_TB ('dotplot', feature_vec = xx [, color_col], 
                        rotation=0, color_fill=T, aes_param=AP) +
                        ggplot2::theme (axis.text.y = element_text (size=amplify_font*AP$fontsize))+
                        ggplot2::guides (fill= ggplot2::guide_legend (override.aes=
                                                        list(color='white', alpha=1)))
        }
        return (plot_xx)
}
