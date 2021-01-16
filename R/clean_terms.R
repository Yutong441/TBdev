# This script should be more appropriately called `ontology.R` because it deals
# with not only KEGG terms, but GO, Reactome, GSEA as well
#
# Install the latest version of clusterProfiler
# It appears that this is not a problem for R 4.0 according to 
# https://github.com/YuLab-SMU/clusterProfiler/issues/245
# However, for R 3.6, I need to:
# remove.packages ('enrichplot')
# remove.packages ('DOSE')
# remove.packages ('fgsea')

# devtools::install_github("ctlab/fgsea")
# devtools::install_github('https://github.com/YuLab-SMU/enrichplot')
# devtools::install_github('https://github.com/YuLab-SMU/DOSE')
# devtools::install_github('https://github.com/YuLab-SMU/clusterProfiler')
# 
# The side effect is that the compareClusterResult object will lack the
# `termsim` attribute, which needs to be calculated manually:
# according to https://rdrr.io/github/GuangchuangYu/enrichplot/man/pairwise_termsim.html:
# library (GOSemSim)
# library (enrichplot)
# library (clusterProfiler)
# kk <- compareCluster(gene_list, fun="enrichKEGG", organism="hsa",
#                     pvalueCutoff=0.05)
# d <- godata('org.Hs.eg.db', ont="BP")
# kk <- pairwise_termsim(kk, method="JC", semData = d)

# A similar issue is reported: https://github.com/YuLab-SMU/clusterProfiler/issues/252

#-------------------------------------------------- 

#' Obtain the species name for KEGG
#'
#' @param name common name for the species
#' @noRd
get_kegg <- function (name){
        org_names <- c('human'='hsa', 'mouse'='mmu', 'marmoset'='cjc')
        if (!name %in% names (org_names)){
                org <- KEGGREST::keggList ('organism')
                kegg_name <- org [grep (name, org [, 'species']), 'organism']
                rm (org)
                return (kegg_name[1])
        }else{
                return (org_names [names (org_names) == name ])
        }
}

remove_terms <- function (xx_res, AP){
        for ( i in AP$remove_keys){
                rm_index <- grepl (i, xx_res, ignore.case=T)
                xx_res <- xx_res [!rm_index]
        }
        return (xx_res)
}

#' Remove terms from GO/KEGG analysis that are not related to development
#'
#' @param xx a compareClusterResult object
#' @noRd
clean_terms <- function (xx, AP, attr_name='compareClusterResult',
                         term_column='Description'){
        xx_res <- attr (xx, attr_name) 
        keep_terms <- remove_terms (xx_res [, term_column], AP)
        xx_res  <- xx_res [xx_res [, term_column] %in% keep_terms, ]
        attr (xx, attr_name) <- xx_res
        return (xx)
}

simplify_terms <- function (xx, simple_terms, term_col='Description',
                            attr_name='compareClusterResult'){
        dat <- attr (xx, attr_name) 
        new_terms <- as.character (simple_terms$sub [match (
                                   dat [, term_col], simple_terms$ori) ])
        na_terms <- is.na (new_terms)
        new_terms [na_terms] <- as.character (dat [na_terms, term_col])
        names (new_terms) <- dat [, term_col]
        new_terms [!na_terms & new_terms %in% dat [, term_col] ] <- NA
        new_terms [!na_terms & duplicated (new_terms) ] <- NA
        new_terms <- gsub ('signaling pathway', 'signaling', new_terms) 
        # Cancer terms are simplified because they are not the focus of this research
        cancer_terms <- grep ('cancer', new_terms, ignore.case=T)
        new_terms [cancer_terms [1] ] <- 'Cancer'
        if (length (cancer_terms)>1) {new_terms [cancer_terms [2:length (cancer_terms)]] <- NA}
        return (new_terms)
}

simplify_gsea <- function (xx_df, simple_terms=NULL, term_col='category'){
        if (!is.null (simple_terms)){
                new_terms <- as.character (simple_terms$sub [match (
                                        xx_df [, term_col], simple_terms$ori) ])
                na_terms <- is.na (new_terms)
                new_terms [na_terms] <- as.character (xx_df [na_terms, term_col])
                names (new_terms) <- xx_df [, term_col]
                new_terms <- gsub ('signaling pathway', 'signaling', new_terms) 
                cancer_terms <- grep ('cancer', new_terms, ignore.case=T)
                new_terms [cancer_terms [1] ] <- 'Cancer'
                xx_df [, term_col] <- new_terms
        }
        return (xx_df)
}

append_default_dictionary <- function (new_dict, append_default){
        if (append_default){
                data (GOsimp)
                if (is.null (new_dict)){new_dict <- GOsimp
                }else{
                        new_dict <- rbind (new_dict, GOsimp) 
                        new_dict <- new_dict [!duplicated (new_dict$ori), ]
                }
        }
        return (new_dict)
}

#' For the results from `compareClusterResult` object
#'
#' @description This function will print all the GO/KEGG/Reactome terms that
#' match the requested keyword or phrases with their p values. It will
#' optionally return a list of all the genes corresponding to the matches terms
#' @param xx a compareClusterResult object
#' @param term a keyword that matches the term of interest
#' @param organism_db which organism database, for example, org.Hs.eg.db
#' @param category_col which column in x that has the terms
#' @param return_val whether to return the gene names
#' @examples
#' KG <- modules::use ('KEGG_path.R')
#' library (GOSemSim)
#' d <- godata('org.Hs.eg.db', ont="BP")
#' kk <- KG$compare_cluster_enrichment (markers, d, enrich_area='KEGG')
#' path_val <- KG$gene_per_term (kk, 'JAK-STAT', return_val=T)
#' @export
gene_per_term <- function (xx, term, organism_db, category_col='Description', return_val=F){
        if (class (xx) == 'compareClusterResult' ){
                xx_df <- xx@compareClusterResult
        }else{ xx_df <- xx }
        all_index <- grep (term, xx_df [, category_col], ignore.case=T)
        all_terms <- list ()
        for (i in 1:length (all_index)){
                print (paste ('the genes for', xx_df[all_index [i], category_col ] ))
                print (paste ('pvalue is', xx_df$p.adjust [all_index [i] ] ) )
                all_genes <- xx_df$geneID [all_index [i] ]
                all_genes <- strsplit (all_genes, '/') [[1]]
                gene_entrez <- AnnotationDbi::mapIds(organism_db, as.character (
                                                all_genes),  'SYMBOL', 'ENTREZID')
                names(gene_entrez) <- NULL
                all_terms[[i]] <- gene_entrez
                print (paste ('there are', length (gene_entrez), 'genes' ) )
        }
        if (return_val){
                names (all_terms) <- xx_df [all_index, category_col]
                return (all_terms)
        }
}

#' Obtain genes of each enriched terms
#'
#' @description For the data frame generated from enrichment analysis, achieves
#' the same function as `gene_per_term`
#' @importFrom magrittr %>%
#' @export
gene_per_enriched_term <- function (xx, term, organism_db, cell_type=NULL, category_col='category'){
        if (!is.null (cell_type) ){
                xx_df <- xx %>% filter (cluster == cell_type)
        }else{ xx_df <- xx }
        all_index <- grep (term, xx_df [, category_col], ignore.case=T)
        uniq_terms <- unique ( xx_df [all_index, category_col] )
        all_terms <- list ()
        for (i in 1:length (uniq_terms) ){
                print (paste ('the genes for', uniq_terms [i] ))
                xx_df %>% dplyr::filter ( !!as.symbol (category_col) %in% uniq_terms[i] ) %>%
                        rownames () -> all_genes
                all_genes <- sapply (all_genes, function (x) {strsplit (x, '\\.')[[1]][3]})

                gene_entrez <- AnnotationDbi::mapIds(organism_db, as.character (
                                                all_genes),  'SYMBOL', 'ENTREZID')
                names(gene_entrez) <- NULL
                all_terms[[i]] <- gene_entrez
        }
        return (all_terms)
}
