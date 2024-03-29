# -------------------------
# Dimensionality reduction 
# -------------------------

#' Run diffusion map
#'
#' @param x Seurat object
#' @param reduction run diffusion map on the reduced dimension of which
#' method
#' @param dims run diffusion map on which reduced dimensions
RunDiffusion <- function (x, reduction='pca', dims=1:10, select_cells=NULL){
        diff_input <- x@reductions[[reduction]]@cell.embeddings [, dims]
        # this function from destiny is 4.23% faster than that in the scater
        # package
        if (is.null (select_cells)){
                dff_out <- destiny::DiffusionMap (as.matrix (diff_input))
                diff_embed <- dff_out@eigenvectors
        }else{
                diff_sel <- diff_input [select_cells, ]
                diff_out_sel <- destiny::DiffusionMap (as.matrix (diff_sel))
                diff_not_sel <- destiny::dm_predict (diff_out_sel, as.matrix (
                                                   diff_input [!select_cells, ]))
                rownames (diff_not_sel) <- rownames (diff_input)[!select_cells]
                diff_embed <- rbind (as.matrix (diff_out_sel@eigenvectors), 
                                     as.matrix (diff_not_sel))
        }
        diff_embed <- diff_embed [match (colnames (x), rownames(diff_embed)), ]
        rownames (diff_embed) <- colnames (x)
        diff_map_seurat <- Seurat::CreateDimReducObject (embeddings= diff_embed, key =
                                                 'DC_', assay='RNA')
        x@reductions[['diff_map']] <- diff_map_seurat
        return (x)
}

#' Add GPLVM data to Seurat dim red assay
#'
#' @param x a Seurat object
#' @param embedding the mean pseudotime estimated from GrandPrix
#' @export
RunGPLVM <- function (x, embedding, assay='RNA'){
        rownames(embedding) <- colnames (x)
        gp_seurat <- Seurat::CreateDimReducObject (embeddings= embedding, 
                                                   key = 'GP_', assay=assay)
        x[['gplvm']] <- gp_seurat
        return (x)
}

#' Run DR
#'
#' @param assay assay to use for dimensionality reduction
#' @param select_cells run PCA on a selected population of cells;
#' should be a boolean vector
#' @param feature_num how many variably expressed genes to compute
#' @param find_var_features whether to run variable gene analysis
#' @param save_mem for large matrices whether the entire matrix should be
#' stored, or simply storing the `VariableFeatures`
#' @param normalize whether to normalize the data
#' @param var_scale the variable features only to reduce memory
#' @param pca_features run PCA on a selected list of genes, should be a
#' character vector
#' @param run_diff_map whether to run diffusion map using the package
#' destiny. Depending on the dataset, this may take several minutes
#' @param run_umap whether to run umap
#' @param cluster whether to perform clustering as well
#' @return a Seurat object
#' @importFrom Seurat VariableFeatures
#' @export
run_dim_red <- function (x, assay=NULL, 
                         select_cells=NULL,
                         feature_num=2000, find_var_features=FALSE, 
                         save_mem=NULL, 
                         normalize=FALSE, 
                         var_scale=FALSE, pca_features=NULL, 
                         run_diff_map=FALSE, run_umap=TRUE, 
                         cluster=TRUE, cluster_res=0.5, neighbor_dim=10){

        if (normalize){x <- Seurat::NormalizeData(object = x, normalization.method =
                         "LogNormalize", scale.factor = 10000)
        }
        # find variable features
        if (is.null(assay)){assay <- Seurat::DefaultAssay (x) }
        Seurat::DefaultAssay (x) <- assay
        var_features <- VariableFeatures (x)
        if (length (var_features) == 0 | is.null (var_features) | find_var_features ){
                x <- Seurat::FindVariableFeatures(x, assay= assay, selection.method =
                                          "vst", nfeatures = feature_num)
        }

        if (!is.null(save_mem)){ 
                orig_var_feature <- VariableFeatures (x)
                x <- Seurat::FindVariableFeatures(x, assay= assay, selection.method =
                                          "vst", nfeatures = save_mem*feature_num)
                x <- x [VariableFeatures (x), ] 
                VariableFeatures (x) <- orig_var_feature
        }
        print_dim (x) # from exp_mat.R
        # dimensionality reduction
        if (var_scale) {x <- Seurat::ScaleData(x, assay=assay, features = VariableFeatures (x))
        }else {x <- Seurat::ScaleData(x, assay=assay, features = rownames(x))}
        x <- selected_cell_PCA (x, select_cells=select_cells, assay=assay,
                                select_genes=pca_features)
        n_umap <- pmin (dim (x)[2]-1, 30 )
        if (run_umap) {x <- Seurat::RunUMAP (x, assay=assay, dims=1:n_umap, n.neighbors=n_umap)}
        if (run_diff_map) {x <- RunDiffusion (x, reduction='pca', dims=1:10, 
                                              select_cells=select_cells)}

        if (cluster){
                x <- Seurat::FindNeighbors(x, dims = 1:neighbor_dim)
                x <- Seurat::FindClusters(x, resolution = cluster_res)
        }
        return (x)
}

#' Run PCA on selected population of cells
#' 
#' @param x Seurat object
#' @param select_cells PCA on a subset of cells, a boolean vector;
#' default is all cells
#' @param select_genes PCA on a subset of genes, default is variable
#' genes
#' @importFrom Seurat VariableFeatures
selected_cell_PCA <- function (x, select_cells=NULL, assay=NULL, select_genes=NULL){
        # initialise default parameters
        if (is.null(assay)){assay <- Seurat::DefaultAssay (x) }
        if (is.null(select_genes)){select_genes <- VariableFeatures(x)}
        num_pc <- pmin ( dim (x)[2]-1, 50 )
        if (!is.null(select_cells)){
                # perform PCA
                select_cells [is.na (select_cells)] <- TRUE
                exp_mat <- t (as.matrix (x[[assay]]@scale.data [select_genes, ]))
                pc_select <- stats::prcomp (exp_mat [as.vector(select_cells),], center=F, scale=F)
                print (dim (pc_select$x))
                pc_nonsel <- exp_mat [!select_cells, ] %*% pc_select$rotation

                # reorganise results
                rownames (pc_select$x) <- rownames (exp_mat [select_cells, ])
                rownames (pc_nonsel) <- rownames (exp_mat [!select_cells, ])
                all_pc <- rbind (pc_select$x, pc_nonsel)
                all_pc <- all_pc [match (rownames (exp_mat), rownames (all_pc)), ]

                # add into seurat object
                pca_seurat <- Seurat::CreateDimReducObject (embeddings= all_pc,
                                                    loadings=pc_select$rotation, key =
                                                            'PC_', assay=assay)
                x@reductions[['pca']] <- pca_seurat
        }else{x <- Seurat::RunPCA (x, assay=assay, features=select_genes, npcs=num_pc)}
        return (x)
}

# ---------- Visualisation ---------- #
label_or_not <- function (x, threshold=9){
        if (length (unique(x)) > threshold){return (FALSE)
        }else{return (TRUE)}
}

label_order <- function (x){
        unique_type <- unique (x)
        return ( unique_type[order (nchar (unique_type), unique_type)])
}

#' Quick way of plotting DimPlot in Seurat
#' 
#' @description Please run the `run_dim_red` function before executing this
#' function. Internally, this function calls `gg_DimPlot` for 2D plotting and
#' `dim_red_3D` for 3D. Please refer to these 2 functions for further parameter
#' settings
#'
#' @param input_data a Seurat object
#' @param group.by a 2-element vector, parsed to the `group.by` argument of
#' `DimPlot` in Seurat
#' @param plot_3D whether to plot in 3D
#' @param DR the type of dimensionality reduction method
#' @param select_cells cells to be shown in the plot
#' @param dims two axes of the reduced dimensional space
#' @param num_col how many columns the subplots are arranged
#' @param no_facet whether to use `ggplot2::facet_wrap`. Currently only 3D
#' version is supported.
#' @return a ggplot or list of ggplot objects
#' @export
plot_dim_red <- function (input_data, group.by, DR='pca', select_cells=NULL,
                          dims=c(1,2), return_sep=F,
                          num_col=NULL, num_row=NULL, no_facet=T, AP=NULL,...){
        if (is.null(select_cells)){
                plot_data <- input_data
        }else{plot_data <- input_data [, select_cells]}

        # determine whether to plot in 2D or 3D
        if (length (dims) == 3 ){plot_3D <- TRUE
        }else{
                plot_3D <- FALSE
                if (length (dims) > 3){
                        dims <- dims [1:2]
                        print ('More than 3 dimensions being supplied, only the
                               first two will be used')
                }
        }
        plots <- list()
        for (i in 1:length(group.by)){
                if (!plot_3D){
                        plots[[i]] <- gg_DimPlot (plot_data, feature=group.by[i], DR=DR, 
                                                     dims=dims, AP=AP, ...)
                }else{
                        if (no_facet){
                                plots [[i]] <- DimPlot_3D (input_data, group.by[i],
                                                           DR=DR, AP=AP,...)
                        }
        }}
        if (!no_facet){
                plots <- DimPlot_3D (input_data, group.by, DR=DR, AP=AP,...)
                return_sep <- T
        }
        if (!return_sep){
                return (ggpubr::ggarrange (plotlist=plots, ncol=num_col, nrow=num_row) )
        }else{return (plots) }
}

#' Plot the contribution of genes to each PCs
#'
#' @param x a Seurat object
#' @param group.by which feature of genes to color, default is the cell
#' type in which the genes are differentially expressed
#' @param dims which pcs to plot
#' @param show_markers the genes to shown in the final plot. If `NULL`, all genes
#' will be shown. If `TF`, a list of transcriptional factors will be shown.
#' Alternatively, the argument accepts a vector of genes.
#' @param directory where the DE gene results are stored. If NULL, the cells
#' will not be coloured.
#' @param label the label of the DE gene data
#' @param cluster_col which column in the DE gene result that stores the
#' cluster names
#' @param feature_col which column in the DE gene result that stores the gene
#' names
#' @param AP aesthetic parameters
#' @importFrom ggplot2 aes aes_string
#' @importFrom magrittr %>%
#' @export
plot_gene_PC <- function (x, group.by=NULL, show_markers=NULL, dims=c(1,2),
                          directory=NULL, label='all', DR='pca',
                          color_markers=NULL, DE_rank='logFC',
                          cluster_col='cluster', feature_col='feature', AP=NULL){
        AP <- return_aes_param (AP)
        plot_PC <- paste ('PC', dims, sep='')
        plot_gene <- x@reductions[[DR]]@feature.loadings 

        # run DE gene analysis
        if (!is.null (directory) & !is.null (group.by)){
                print ('finding marker genes')
                markers <- find_DE_genes (x, directory, group.by=group.by, label=label)
                markers %>% dplyr::arrange (desc (!!as.symbol (DE_rank))) -> markers
                markers <- markers [match (rownames (plot_gene), markers [, feature_col]), ]
        }else{markers <- data.frame (cluster=rep ('all', nrow(plot_gene) ))}

        print ('select markers')
        if (is.null (show_markers)){marker_list <- rownames (plot_gene)
        }else if (show_markers=='TF'){
                #data (TF, package='TBdev')
                data (TF)
                marker_list <- TF
        }else{marker_list <- show_markers}

        markers [, cluster_col ] <- partial_relevel (markers [, cluster_col],
                                                     AP$cell_order)

        if (is.null (color_markers)){  color_cluster <- markers [, cluster_col]
        }else{ color_cluster <- names (color_markers) [  match (
                                        rownames (plot_gene), color_markers) ] }

        print ('plotting')
        plot_gene %>% as.data.frame () %>%
                magrittr::set_colnames (paste ('PC', 1:dim(plot_gene)[2], sep='')) %>%
                dplyr::select (plot_PC) %>%
                tibble::add_column (gene = rownames (plot_gene)) %>%
                tibble::add_column (cluster = color_cluster) %>%
                dplyr::filter (gene %in% marker_list) -> plot_data
        if (!is.null(group.by)){
                plot_data %>% dplyr::filter (cluster %in% unique (x@meta.data[, group.by])) -> plot_data
        }

        ggplot2::ggplot (plot_data, aes_string(plot_PC[1], plot_PC[2]))+
                ggplot2::geom_text (aes(label=gene, color=cluster), size=AP$point_fontsize, 
                                    family=AP$font_fam) -> plot_ob
        return (plot_ob +theme_TB ('dim_red', plot_ob=plot_ob, feature_vec=plot_data$cluster, 
                                   aes_param=AP))
}
