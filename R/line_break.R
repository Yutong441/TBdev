order_genes_one_cond <- function (sel_markers, gdf){
        sel_markers [match (gdf$gene, sel_markers$feature),] %>%
                dplyr::select (logFC) %>% cbind (gdf) %>%
                dplyr::arrange (dplyr::desc (abs (logFC)))
}

order_genes <- function (genes_df, markers, compare_name){
        clusters <- unique (genes_df$celltype)
        markers <- markers [markers$compare_group == unique (compare_name),]
        if (mean (clusters %in% markers$group) == 1){
                glist <- lapply (as.list(clusters), function(x){
                        order_genes_one_cond (markers[markers$group==x,],
                                              genes_df[genes_df$celltype==x,]
                        )})
                return (do.call (rbind, glist))
        }else{
                print ('the DE gene marker dataframe is not obtained from the same
                       clusters as the GSEA.')
                return (genes_df)
        }
}

#' Obtain the label information for enrichment barplot
#'
#' @description This function joins the enrichment term with the top
#' differentially expressed genes in a single string, that will be shown in the
#' enrichment barplot labels. 
#' @param xx raw dataframe from GSEA
#' @param sum_gsea a dataframe generated from `summarise_gsea`
#' @param markers a dataframe generated from `find_DE_genes`. If NULL, the
#' genes will be ordered alphabetically.
#' @return a character vector of labels
#' @importFrom magrittr %>%
term_gene_labels <- function (xx, sum_gsea, organism_db, show_gene_labels=3,
                              markers=NULL, category_col= 'category',
                              cluster_col='cluster', term_col='termID'){
        cat_col <- sum_gsea [, category_col]
        clust_col <- sum_gsea [, cluster_col]
        text_df <- xx [xx [, category_col] %in% cat_col, ]

        term_genes <- lapply (as.list(1:nrow(sum_gsea)), function (i){ 
                text_df [text_df [, category_col] == cat_col[i] & text_df [
                          , cluster_col] == clust_col[i], term_col]  })
        cat_clust <- paste (cat_col, '___', clust_col, sep='')
        names (term_genes) <- paste (cat_clust, '__', sep='')

        # merge all entrez ID into a single vector, which is converted into
        # common gene names. This is much faster.
        genes_vec <- gsea_entrez_to_name (do.call (c, term_genes), organism_db)
        names (genes_vec) <- gsub ('__[0-9]*$', '', names (genes_vec) )
        genes_df <- gene_vec_to_df (genes_vec)

        if (!is.null (markers) ){
                compare_name <- unique (xx$compare_group)
                genes_df <- order_genes (genes_df, markers, compare_name)
        }
        term_genes <- lapply (as.list (cat_clust), function (x) {paste (line_break_every(
                                genes_df[genes_df$meta==x, 'gene'][1:show_gene_labels]), 
                                collapse=', ')}) %>% unlist()
        # remove the NA terms that would appear if the number of enriched genes
        # is smaller than `show_gene_labels`
        term_genes <- gsub (', NA', '', term_genes)
        # prevent extra commas from appearing after each line break
        term_genes <- gsub ('\n,', '\n', term_genes)
        # prevent empty spaces
        term_genes <- gsub (',( \n)+$', '', term_genes)

        return (paste (cat_col, '\n(', term_genes, ')', sep=''))
}

line_break_every <- function (vec, separator='\n', every_n=2){
        vec <- as.character (vec)
        n_interval <- floor(length(vec)/every_n)
        total_n <- n_interval*(every_n+1) + length(vec)%%every_n
        final_vec <- rep (separator, total_n)
        final_ind <- rep (T, total_n)
        final_ind [seq (every_n+1, total_n, by=every_n+1)] <- F

        final_vec [final_ind] <- vec
        return (final_vec)
}

strsplit_add_sep <- function (string, sep){
        vec <- str_split (string, sep)[[1]]
        if (length(vec)>1){vec <- paste (vec, sep,sep='')}
        vec[length(vec)] <- gsub (sep, '', vec[length(vec)])
        return (vec)
}

#' @importFrom magrittr %>%
wrap_one_long_sentence <- function (sentence, index_sepa, break_sepa, thres){
        # break up a sentence into lines that already exist
        str_x <- strsplit (sentence, break_sepa)[[1]]
        sapply (str_x, function (i){
                # if the sentence is long enough
                if (length (strsplit(i,'')[[1]]) > thres){
                        # for each line, search for points that 
                        # allow for breaking, e.g. commas
                        str_i <- strsplit_add_sep (i, index_sepa)
                        sapply (1:length(str_i), function (k) {
                                all_len <- length (strsplit (str_i[k], '')[[1]])
                                # need to break if the next word is too long
                                if (k + 1 <= length (str_i)) {
                                        all_len2 <- length (strsplit (str_i[k+1], '')[[1]])
                                }else{all_len2 <- all_len}
                                if (all_len > thres | all_len2 > thres){
                                        paste (str_i[k], break_sepa, sep='')
                                }else{str_i[k]}
                        # because `strsplit_add_sep` above, there is no need to
                        # collapse with a separator
                        }) %>% paste (collapse='') %>% return ()
                # if the sentence is short, there is no need for breaking
                }else{return (i)}
        }) %>% paste (collapse=break_sepa) %>% return ()
}

#' Line break after very long words
#'
#' @param vec a character vector, e.g. `c('HLA-G','TBX3','CGA')`
#' @param index_sepa where may a line break be inserted
#' @param break_sepa the separator for line break, i.e. '\n' for most
#' applications except `<br>` for markdown syntax.
#' @param thres how many letters maximum in each line
#' @description I could have used `stringr::str_wrap` but that function does
#' not consider if the input string already has breaks
wrap_long_sentence <- function (vec, index_sepa=',', break_sepa='\n', thres=8){
        sapply (vec, function(x){
                wrap_one_long_sentence (x, index_sepa, break_sepa, thres)
        })
}

line_break <- function (vec, separator='\n', every_n=2, thres=8){
        vec <- paste (line_break_every(vec, separator, every_n), collapse=',')
        vec <- wrap_long_sentence (vec, break_sepa=separator, thres=thres)
        # remove duplicated line breaks, i.e. don't want any empty lines
        vec <- gsub ('\n\n','\n', vec)
        # remove duplicated commas
        vec <- gsub (',,',',', vec)

        # remove NA fields
        vec <- gsub (', NA', '', vec)
        # prevent extra commas from appearing after each line break
        vec <- gsub ('\n,', '\n', vec)
        # prevent empty spaces
        vec<- gsub (',( \n)+$', '', vec)
        return (vec)
}

deparse_labels <- function (f1, f2){
        deparse (bquote (atop (bold (.(f1)), .(f2))))
}
