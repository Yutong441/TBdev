# follow the method from Liu 2018, which follows Tirosh 2016
# Tirosh 2016 follows the method of Masosko 2015
# download the genes from Macosko 2015 Table S2
# references:
# Liu 2018: https://www.nature.com/articles/s41422-018-0066-y
# Tirosh 2016: https://science.sciencemag.org/content/352/6282/189  
# Macosko 2015: https://www.sciencedirect.com/science/article/pii/S0092867415005498#mmc2
# Method of Macosko 2015: https://jdblischak.github.io/singleCellSeq/analysis/cell-cycle.html

# -------------------
# cell cycle analysis
# -------------------
get_cycle_order <- function(){
        return (c('M', 'M.G1', 'G1.S', 'S', 'G2.M'))
}
load_cell_cycle_genes <- function (){
        data (cell_cycle)
        gene_list <- list ()
        for (i in 1:ncol(cell_cycle)){
                gene_list [[i]] <- as.vector (unique ( cell_cycle[, i] ))
                names (gene_list [[i]] ) <- rep ( colnames (cell_cycle)[i], 
                                                 length (gene_list[[i]]) )
        }
        gene_vec <- do.call ( c, gene_list )
        gene_vec <- trimws (gene_vec, which='both') #remove white space
        return (gene_vec [gene_vec != '']) #remove empty fields
}

#' Assign phase score
#'
#' @description To calculate the phase score, Macosko (2015) used mean
#' expression of genes that mark a particular phase of the cell cycle. Not all
#' genes in the list are used. Only those that correlate with the mean
#' expression are included, followed by some normalisation. I am not certain
#' whether to apply normalisation or not, given that this research uses batch
#' combined data that has already been normalised. Qualitatively in a heatmap,
#' the results did not look much different.
#' @param x an expression matrix
#' @param gene_vec a vector of genes with their names being which cell cycles
#' they are expressed
#'
#' @export
get_phase_score <- function (x, gene_vec = NULL, norm_factor=F, norm_column=F, norm_row=F){
        if (is.null (gene_vec)){gene_vec <- load_cell_cycle_genes ()}
        # convert gene vector into list to use lapply
        gene_list <- lapply ( as.list (unique (names (gene_vec)) ), 
                             function (x) {gene_vec [ names (gene_vec) %in% x ]} )
        names (gene_list) <- unique ( names (gene_vec) )

        ans <- lapply(gene_list, function(xx){
                reads_single_phase <- x[rownames(x) %in% unlist(xx) ,]
                # add average expression
                combined_matrix <- rbind(reads_single_phase,average=apply(reads_single_phase,2,mean))
                cor_matrix <- cor(t(combined_matrix))
                cor_vector <- cor_matrix[,dim(cor_matrix)[1]]

                # restrict to correlation >= 0.3 
                reads_single_phase_restricted <- reads_single_phase[rownames(reads_single_phase) 
                                                %in% names(cor_vector[cor_vector >= 0.3]),]
                # remove columns with all zeros
                filter_out <- colMeans (reads_single_phase_restricted) == 0.
                reads_single_phase_restricted <- reads_single_phase_restricted [,!filter_out]

                # apply normalization to reads
                if (norm_factor){
                        norm_factors_single <-
                                edgeR::calcNormFactors(reads_single_phase_restricted,
                                                method = "TMM")
                        lib_size <- colSums(x)[!filter_out] * norm_factors_single
                        reads_single_cpm <- SingleCellExperiment::cpm(
                        reads_single_phase_restricted, log = T, lib.size = lib_size)
                        #output the phase specific scores (mean of normalized
                        #expression levels in the phase)
                }else{
                        reads_single_cpm <- reads_single_phase_restricted 
                }
                return (apply(reads_single_cpm,2,mean))
        })
        cell_names <- lapply (ans, function(x){names (x)})
        shared_cells <- Reduce (intersect, cell_names)
        ans <- lapply (ans, function (x) { x [match (shared_cells, names (x) )] } )
        ans <- do.call (cbind, ans)

        if (norm_column){ans <- scale (ans) }
        if (norm_row){ ans <- t( scale ( t( ans ) ) ) }
        return (ans)
}

#' Make a heatmap for phase score
#'
#' @param x a matrix containing the phase scores
#' @param seurat_ob a seurat object that provides the metadata
#' @param group.by which feature to plot the row color bar for the heatmap
#' @export
make_cycle_heat <- function (x, seurat_ob, group.by=NULL, AP=NULL, ...){
        cycle_names <- factor ( colnames (x), levels=get_cycle_order() )
        print (cycle_names)
        x_seurat <- Seurat::CreateSeuratObject (x, meta.data=data.frame (
                        cell_cycle=cycle_names, row.names=colnames (x)) )
        cell_vec <- rownames (x_seurat)
        if (!is.null (group.by)){
                names (cell_vec) <- seurat_ob@meta.data [match (rownames (x), 
                                                colnames (seurat_ob)), group.by]
        }
        seurat_heat (x_seurat, color_row=cell_vec, group.by='cell_cycle',
                     slot_data='counts', show_row_names=F, cluster_columns=F,
                     row_scale=F, cluster_rows=T, heat_name = 'phase score', 
                     column_legend_labels = 'cell cycle',
                     row_legend_labels = 'cell type', center_scale=T,
                     column_reorder_levels=list (get_cycle_order()), AP=AP, ...)
}

#' Correlation of phase score with PC
#' 
#' @param x a seurat object
#' @param num_dim which dimension in the dimensionality reduction (DR) method to use
#' @param corr_feature which feature is the selected DR correlated with
#' @importFrom ggplot2 aes aes_string
#' @examples
#' exp_mat <- as.matrix (Seurat::GetAssayData (all_data, assay='RNA', slot='data'))
#' ans <- get_phase_score (exp_mat )
#' phase_PC_plot (ans, all_data, 'broad_type')
#' @export
phase_PC_plot <- function (ans, x, corr_feature, DR='pca', num_dim=1, AP=NULL){
        AP <- return_aes_param (AP)
        pca <- x@reductions[[DR]]@cell.embeddings 
        pca <- pca [match (rownames (ans), rownames (pca)), ]

        metadata <- x@meta.data [match (rownames (ans), colnames (x) ), ]
        data.frame (ans) %>%
                tibble::add_column (Type=metadata [, corr_feature]) %>%
                tibble::add_column (PC1=pca [, num_dim]) %>%
                tidyr::gather ('cycle_phase','phase_score', -Type, -PC1) %>%
                dplyr::mutate (cycle_phase = factor (cycle_phase, levels=c(
                                        'M', 'M.G1', 'G1.S', 'S', 'G2.M') ) ) -> plot_data
        ggplot2::ggplot (plot_data, aes_string (x='PC1', y='phase_score')) +
                ggplot2::geom_point (aes (fill=Type),
                                     color=AP$point_edge_color,
                                     shape=AP$normal_shape,
                                     size=AP$pointsize) +
                ggplot2::geom_smooth (method='lm') +
                ggpubr::stat_cor(aes(label= ..rr.label..),
                                 label.x.npc='left',
                                 size=AP$point_fontsize) +
                ggplot2::xlab (colnames (pca)[num_dim] )+
                ggplot2::facet_wrap (~cycle_phase)+ ggplot2::labs (fill='')+
                theme_TB ('dotplot', feature_vec=metadata [, corr_feature], 
                          aes_param=AP, rotation=0, color_fill=T)+
                custom_tick (plot_data$phase_score)

}

#' Plot the distribution of cycle phase scores as ridge/wave like plots
#'
#' @description same usage as `phase_PC_plot`
#' @importFrom ggplot2 guide_legend aes
#' @export
phase_density_plot <- function (ans, x, corr_feature, select_phase=NULL,
                                num_row=2, AP=NULL){
        AP <- return_aes_param (AP)
        metadata <- x@meta.data [match (rownames (ans), colnames (x) ), ]
        if (is.null (select_phase)){select_phase <- colnames (ans) }
        data.frame (ans) %>%
                tibble::add_column (Type=metadata [, corr_feature]) %>%
                tidyr::gather ('cycle_phase','phase_score', -Type) %>%
                dplyr::filter (cycle_phase %in% select_phase) %>%
                dplyr::mutate (cycle_phase = factor (cycle_phase, levels = 
                                        get_cycle_order()) ) -> plot_data
        ggplot2::ggplot (plot_data, aes (x=phase_score, y=..density.., fill=Type)) +
                ggplot2::geom_density (alpha=AP$ridge_alpha)+
                ggplot2::theme_minimal () + ggplot2::labs (fill='')+
                ggplot2::facet_wrap (~cycle_phase, nrow=num_row) + 
                theme_TB ('dotplot', feature_vec = metadata[, corr_feature],
                          color_fill=T, rotation=0, AP=AP)+
                ggplot2::guides(fill = guide_legend(override.aes = list(alpha = 1, 
                                                        color=AP$point_edge_color))) 
}

#' Plot cycle score over pseudotime
#' 
#' @param plot_data a dataframe with cell cycle scores and pseudotime
#' information
#' @param time_col where the pseudotime information is stored
#' @param sel_phase which cell cycle phases to display, can be M, M.G1, G1.S, S, G2.M
#' @param color_by which column to color
#' @param num_col in how many columns the subplots are arranged
#' @importFrom ggplot2 aes aes_string
#' @importFrom magrittr %>%
#' @examples
#' exp_mat_cc <- as.matrix (GetAssayData (all_data, assay='RNA', slot='data'))
#' ans <- get_phase_score (exp_mat_cc)
#' plot_data <- cbind (ans, all_data@meta.data)
#' cycle_over_time (plot_data, time_col='MGP_PT')
#' @export
cycle_over_time <- function (plot_data, color_by, time_col,sel_phase=c('M', 'S'), 
                             num_col=1, extend_max=0.4, AP=NULL){
        AP <- return_aes_param (AP)
        plot_data %>% dplyr::select ( all_of ( c ( sel_phase, color_by, time_col )  ) ) %>%
                tidyr::gather ('phase', 'score', -!!as.symbol (color_by), 
                               -!!as.symbol (time_col)) -> plot_data
        ggplot2::ggplot (plot_data, aes_string (x=time_col, y= 'score')) + 
                ggplot2::geom_point (aes_string (fill=color_by), shape=21, size=AP$pointsize, 
                                     color=AP$point_edge_color) +
                ggpubr::stat_cor(aes (label= ..rr.label..), label.x.npc='center', 
                                 size=AP$point_fontsize) +
                ggplot2::facet_wrap (~phase, ncol=num_col) +ggplot2::labs (fill='')+
                custom_tick (plot_data$score) + 
                custom_scale (plot_data [, time_col], 'x', extend_ratio_max=extend_max) +
                theme_TB ('dotplot', feature_vec = plot_data [, color_by],
                          color_fill=T, rotation=0, aes_param=AP)
}
