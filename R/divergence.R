# ----------KL divergence----------
KL_normal <- function (mu1, mu2, var1, var2){
    KL = -0.5 * log(var1) + 0.5 * log(var2) - 0.5 
    KL = KL+ 0.5 * ((mu1 - mu2)^2 + var1) / var2 
    return (KL)
}

WD_normal <- function (mu1, mu2, var1, var2){
        WD <- (mu1 - mu2)^2 + var1 + var2 - sqrt (var1*var2)
        return (WD)
}

sep_mean_val <- function (x){
        all_genes <- colnames (x) [ grep ('^mean_', colnames (x) ) ]
        all_genes <- gsub ('^mean_', '', all_genes)
        pred_mean <- x[, paste ('mean_', all_genes, sep='') ]
        pred_var <- x [, paste ('var_', all_genes, sep='') ]
        colnames (pred_mean) <- all_genes
        colnames (pred_var) <- all_genes
        rownames (pred_mean) <- rownames (x)
        rownames (pred_var) <- rownames (x)
        return (list(as.matrix (pred_mean), as.matrix (pred_var) ))
}

DE_by_KL <- function (pred_list, index1, index2, divergence='KL'){
        if (divergence=='KL'){
                all_KL <- KL_normal (pred_list[[1]][index1, ], pred_list[[1]][index2, ],
                                     pred_list[[2]][index1, ], pred_list[[2]][index2, ])
        }
        if (divergence == 'WD'){
                all_KL <- WD_normal (pred_list[[1]][index1, ], pred_list[[1]][index2, ],
                                     pred_list[[2]][index1, ], pred_list[[2]][index2, ])
        }
        return (all_KL)
}

#' Differential branching gene trajectoy analysis 
#'
#' @param x a dataframe produced by my BRGP script, containing the columns:
#' 'branch' for branch information, 'x' for pseudotime, 'mean_' and 'var_' for
#' the mean and variance of each gene
#' @param num1 the name of one of the branches to compare
#' @param num2 similar to `num1`
#' @param divergence which divergence method to quantify the difference between
#' the 2 distributions. Only support 'WD' (2-Wasserstein distance) and 'KL'
#' (Kullback-Leibler divergence)
#' @importFrom magrittr %>%
#' @export
get_DE_from_KL <- function (x, num1, num2, divergence='WD'){
        rownames (x) <- paste ('cell', 1:nrow(x), sep='')
        pred_list <- sep_mean_val (x)
        index1 <- x['branch'] == num1
        index2 <- x['branch'] == num2

        # compare average fold change
        pred_list[[1]] [ index1 | index2, ] %>% t() -> exp_mat
        metadata <- data.frame ('branch'= x$branch [index1 | index2] )
        rownames (metadata) <- colnames (exp_mat)
        pred_seurat <- Seurat::CreateSeuratObject (exp_mat, meta.data=metadata)
        markers <- presto::wilcoxauc (pred_seurat, 'branch', assay='counts')

        # obtain divergence
        DE_df <- DE_by_KL (pred_list, index1, index2, divergence)
        DE_rank <- colSums (DE_df)
        DE_rank %>% data.frame () %>% magrittr::set_colnames ('divergence') %>% 
                tibble::add_column (gene = names (DE_rank) ) %>% 
                dplyr::arrange (dplyr::desc (divergence) ) -> DE_rank

        # merge information
        markers <- markers [markers$group == num1, ]
        markers$divergence <- DE_rank$divergence [match (markers$feature, DE_rank$gene) ]
        return (markers)
}
