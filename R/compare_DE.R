# ==========Compare between 2 DE lists ==========
# ----------Genes that intersect between the 2 lists----------

get_2_by_2 <- function (vec_true, vec_test, universe){
        out_mat <- matrix (NA, nrow=2, ncol=2)
        rownames (out_mat) <- c('pos', 'neg')
        colnames (out_mat) <- c('pos', 'neg')
        out_mat [1,1] <- length (intersect (vec_true, vec_test))
        out_mat [2,2] <- length (universe [!universe %in% vec_true & !universe %in% vec_test])
        out_mat [1,2] <- length (vec_test [!vec_test %in% vec_true])
        out_mat [2,1] <- length (vec_true [!vec_true %in% vec_test])
        return (out_mat)
}

vivo_vitro_DE <- function (vitro_mark, vivo_mark, vitro1, vitro2, vivo1, vivo2,
                           return_genes=TRUE){
        vitro_mark %>% dplyr::filter (group==vitro1 & compare_group==vitro2 & pval < 0.05 & logFC > 0.5) %>%
                dplyr::pull (feature) -> sel_gene

        vivo_mark %>% dplyr::filter (group==vivo1 & compare_group==vivo2 & pval < 0.05 & logFC > 0.5) %>%
                dplyr::pull (feature) -> sel_gene_vivo

        if (return_genes){
                sect_genes <- intersect (sel_gene, sel_gene_vivo)
                vitro_mark %>% dplyr::filter (feature %in% sect_genes) %>%
                        dplyr::filter (group==vitro1 & compare_group==vitro2) %>%
                        dplyr::arrange (dplyr::desc (logFC) ) %>%
                        dplyr::pull (feature) %>% return ()
        }else{
                get_2_by_2 (sel_gene, sel_gene_vivo, rownames (all_data)) %>% return ()
        }
}

one_vivo_all_vitro <- function (vitro_mark, vivo_mark, vitro_list, vivo1, vivo2){
        N <- length (vitro_list)
        gene_table <- data.frame (vitro1=NA, vitro2=NA, vivo1=NA, vivo2=NA,
                                  TP=NA, FP=NA, FN=NA, gene=NA)
        entry <- 0
        for (i in 1:N){
                for (j in 1:N){
                        if (vitro_list[i] != vitro_list [j]){
                                entry <- entry + 1
                                gene_table [entry, 'vitro1'] <- vitro_list [i]
                                gene_table [entry, 'vitro2'] <- vitro_list [j]
                                vivo_vitro_DE (vitro_mark, vivo_mark,
                                               vitro_list[i], vitro_list[j],
                                               vivo1, vivo2
                                               ) %>% paste (collapse=', '
                                               ) -> gene_table [entry, 'gene']

                                vivo_vitro_DE (vitro_mark, vivo_mark,
                                               vitro_list[i], vitro_list[j],
                                               vivo1, vivo2, 
                                               return_genes=F) -> stat_comp
                                gene_table [entry, 'TP'] <- stat_comp [1,1]
                                gene_table [entry, 'FP'] <- stat_comp [2,1]
                                gene_table [entry, 'FN'] <- stat_comp [1,2]
                        }
                }
        }
        gene_table$vivo1 <- vivo1
        gene_table$vivo2 <- vivo2
        return (gene_table)
}

#' Compare 2 DE lists between all possible combination of cells
#' 
#' @param vitro_mark DE gene table between all possible pairwise comparisons in
#' the in vitro cell types
#' @param vivo_mark DE gene table between all possible pairwise comparisons in
#' the in vivo cell types
#' @param vitro_list all in vitro cells to be compared with each other
#' @param vivo_list all in vivo cells to be compared with each other
#' @export
all_vivo_all_vitro <- function(vitro_mark, vivo_mark, vitro_list, vivo_list){
        N <- length (vivo_list)
        entry <- 0
        mat_list <- list ()
        for (i in 1:N){
                for (j in 1:N){
                        if (i!=j){
                                entry <- entry + 1
                                one_vivo_all_vitro (vitro_mark, vivo_mark,
                                                    vitro_list, vivo_list[i],
                                                    vivo_list[j]) -> mat_list [[entry]]
                        }
                }
        }
        return (do.call (rbind, mat_list))
}

# ----------p-values of the intersecting genes----------
p_val_one_pair <- function (comp_list, vitro_pair, vivo_pair){
        comp_list %>% dplyr::filter (vivo1 %in% vivo_pair & vivo2 %in% vivo_pair &
                                     vitro1 %in% vitro_pair & vitro2 %in% vitro_pair) %>%
        dplyr::select (vivo1, vivo2, vitro1, vitro2, TP) %>% 
        tidyr::unite ('vivo', c('vivo1', 'vivo2'), sep='_vs_') %>%
        tidyr::unite ('vitro', c('vitro1', 'vitro2'), sep='_vs_') %>%
        tidyr::spread (vitro, TP) %>% 
        tibble::column_to_rownames ('vivo') %>%
        as.matrix () %>% return ()
}

gen_pair_list <- function (pair_list, orderless=F){
        N <- length (pair_list)
        if (!orderless){
                pair_mat <- matrix (NA, nrow=N, ncol=N)
                pair_mat [lower.tri (pair_mat, diag=F)] <- 1
        }else{
                pair_mat <- matrix (1, nrow=N, ncol=N)
                diag (pair_mat) <- NA
        }
        pair_vec <- list()
        entry <- 0
        for (i in 1:N){
                for (j in 1:N){
                        if (!is.na (pair_mat [i, j])){
                                entry <- entry + 1
                                pair_vec [[entry]] <- c(pair_list[i], pair_list[j])
                        }
                }
        }
        return (pair_vec)
}

#' Comparison p-value between 2 groups of DE genes
#'
#' @param comp_list output from `all_vivo_all_vitro `
#' @param vitro_list all in vitro cells to be compared with each other
#' @param vivo_list all in vivo cells to be compared with each other
#' 
#' @description
#' First, select the true positive (TP) of the intersection between the following
#' permutations:
#' 
#' | TP       | vitro1 up | vitro2 up |
#' |----------|-----------|-----------|
#' | vivo1 up | a         | b         |
#' | vivo2 up | c         | d         |
#' 
#' The null hypothesis is that the upregulated genes in vitro1 compared with
#' vitro2 (denoted as vitro1_vs_vitro2) do not intersect with vivo1_vs_vivo2 more
#' than vivo2_vs_vivo1, compared with vitro2_vs_vitro1.  
#' 
#' For example, if vitro1 and vitro2 are very similar, then they should both
#' contain a similar number of vivo1_vs_vivo2, and/or vivo2_vs_vivo1. If vitro1
#' and vitro2 are different in the direction of vivo1 and vivo2, then vitro1
#' should contain a different vivo1 to vivo2 ratio, compared with vitro2.
#' 
#' The probability of evaluating this null hypothesis is evaluated by Fisher's
#' exact test:
#' 
#' $$ p = \frac {(a+b)! (c+d)! (a+c)! (b+d)! }{a! b! c! d! (a+b+c+d)!} $$
#'
#' @return a dataframe with the following attributes:
#' `model`: the cell types in comparison. For example, 'bTSC_vs_hTSC-OKAE,
#' ICM_vs_TB' assumes that bTSC is similar to ICM/TB while hTSC-OKAE is similar to
#' the opposite TB/ICM.
#' `ratio`: a/c (see the Method section below)
#' `rev_ratio`: b/d
#' `pval`: Fisher's exact test p value
#' `padj`: adjusted p value for false discovery rate
#' @export
p_val_all_pairs <- function (comp_list, vitro_list, vivo_list){
        vitro_pairs <- gen_pair_list (vitro_list)
        vivo_pairs <- gen_pair_list (vivo_list)
        p_val_df <- data.frame (model=NA, ratio=NA, rev_ratio=NA, pval=NA)
        entry <- 0
        for (i in vitro_pairs){
                for (j in vivo_pairs){
                        entry <- entry +1
                        TP_mat <- p_val_one_pair (comp_list, i, j)
                        TP_ratio <- TP_mat [1,]/TP_mat [2,]
                        p_val_df[entry, 'model'] <- paste (colnames (TP_mat)[1], 
                                                         rownames (TP_mat)[1], sep=', ')
                        p_val_df [entry, 'ratio'] <- TP_ratio [1]
                        p_val_df [entry, 'rev_ratio'] <- TP_ratio [2]
                        p_val_df [entry, 'pval'] <- stats::fisher.test (TP_mat)$p
                }
        }
        p_val_df$padj <- stats::p.adjust (p_val_df$pval, method='fdr')
        return (p_val_df)
}

#' Filter significant pairs of DE genes
#'
#' @param p_val_df output from `p_val_all_pairs `
#' @return a matrix-like dataframe
#' Column names: upregulated genes in vitro1 compared with vitro2
#" Row names: upregulated genes in vivo1 compared with vivo2
#" 
#" forward: vitro1 and vivo1 are similar, vitro2 and vivo2 are similar
#" reverse: vitro1 and vivo2 are similar, vitro2 and vivo1 are similar
#" insignificant: no significant relationship in found
#' @export
signif_p_val <- function (p_val_df){
        all_vivo <- unique (gsub ('^.*,', '', p_val_df$model))
        all_vitro <- unique (gsub (',.*$', '', p_val_df$model))
        p_val_df %>% dplyr::filter (pval < 0.05) %>%
                dplyr::filter (ratio > 1.5 | rev_ratio > 1.5 ) %>%
                dplyr::filter (rev_ratio < 1/1.5 | ratio < 1/1.5) -> filtered_p
        all_vitro %>% as.list () %>% lapply (function (vitro_pair){
                filtered_p %>% dplyr::filter (grepl (vitro_pair, model)) %>%
                        dplyr::mutate (property = ifelse (ratio > rev_ratio, 'forward', 'reverse')) %>%
                        dplyr::mutate (vivo = gsub ('^.*,', '', model) ) %>%
                        dplyr::select (vivo, property) 
        }) -> sum_table_list
        sum_table <- sum_table_list [[1]]
        colnames (sum_table)[2] <- all_vitro [1]
        for (i in 2:length (sum_table_list)){
                sum_table <- merge (sum_table, sum_table_list[[i]], all=T, by='vivo') 
                colnames (sum_table)[i+1] <- all_vitro [i]
        }
        sum_table [is.na (sum_table)] <- 'insignificant'
        return (sum_table)
}

# ----------Pure in vitro DE genes----------
sel_signif_gene <- function (markers, gene_list, comparison, top_num=20){
        if (is.na (top_num)){return (gene_list)
        }else{
                markers %>% dplyr::filter (group == comparison[1] & 
                                           compare_group== comparison[2]) %>%
                        dplyr::filter (feature %in% gene_list) %>%
                        dplyr::slice_max (logFC, n=top_num) %>%
                        dplyr::pull (feature)
        }
}

pure_vitro_DE_genes <- function (vitro_mark, vivo_mark, vitro_list, vivo_list, select_num=NA){
        vitro_pairs <- gen_pair_list (vitro_list, orderless=T)
        vivo_pairs <- gen_pair_list (vivo_list, orderless=T)

        vitro_pairs %>% lapply (function (xy) {
                vitro_mark %>% dplyr::filter (group==xy [1] & compare_group==xy [2] & 
                                              pval < 0.05 & logFC > 0.5) %>%
                        dplyr::pull (feature) -> sel_gene
                vivo_pairs %>% lapply (function (x){
                        vivo_mark %>% dplyr::filter (group==x[1] & compare_group==x[2] & 
                                                     pval < 0.05 & logFC > 0.5) %>%
                                dplyr::pull (feature) }
                ) -> all_vivo_DE
                do.call (c, all_vivo_DE) %>% unique () -> all_vivo_DE
                return (list (length (sel_gene), sel_gene [!sel_gene %in% all_vivo_DE]))
        }) -> pure_vitro
        pure_vitro %>% sapply (function (x){x[[1]]}) -> ori_vitro_num
        pure_vitro %>% lapply (function (x){x[[2]]}) -> pure_vitro
        vitro_pairs %>% sapply (function (x){paste (x, collapse='_vs_')}) -> names (pure_vitro)
        1:length (vitro_pairs) %>% as.list () %>% sapply (function (num_ID){
                sel_signif_gene (vitro_mark, pure_vitro[[num_ID]], 
                                 vitro_pairs [[num_ID]], top_num=select_num) %>% 
                        paste ( collapse=', ')}
        ) -> pure_vitro_collapse
        data.frame (comparison= names (pure_vitro),
                    num=sapply (pure_vitro, length), ori_num= ori_vitro_num,
                    gene=pure_vitro_collapse) -> gene_df
        gene_df %>% dplyr::mutate (retained= num/ori_num) %>% 
                dplyr::relocate (retained, .after=ori_num) %>% return ()
}
