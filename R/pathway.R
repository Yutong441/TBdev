# pathway analysis using module scores

#' Obtain the gene in a particular pathway
#'
#' @param pathway a dataframe with `pathway` column being the comon pathway
#' names, and the `kegg_id` column being the KEGG ID for the pathways
#' @param path_name the common name of a pathway
get_path_genes <- function (pathway, path_name) {
        if (path_name %in% pathway$pathway){
                path_genes <- KEGGREST::keggGet(pathway [ pathway$pathway == path_name, 'kegg_id'])[[1]]$GENE
                print (paste ('found', length (path_genes), 'genes'))
                if (length (path_genes) != 0 ){
                        path_genes <- path_genes [seq ( 2, length (path_genes), 2 ) ]
                        one_path <- sapply (path_genes, function (x) {strsplit (x, ';') [[1]][1]  })
                        return (as.character (one_path))
                }
        }else{
                print ('queried pathway not in the compiled list')
                print ('choose one from below: ')
                print (paste (pathway$pathway, collapse=', ') )
        }
}

heat_high_not <- function (x, color_row, highlight_cells, group.by=c('revised', 'date')){
        if ( is.null (highlight_cells) ){
                print (seurat_heat (x, color_row=color_row, group.by = group.by, 
                                        slot='data', column_rotation=90))
        }else{
                print (seurat_heat_highlight (x, highlight_cells, color_row, group.by, average=T))
        }
}

#' Save the PCA and heatmap of genes in a pathway
#'
#' @description The results will be stored in a pdf file. For each pathway, a
#' PC plot using the genes in that pathway is shown, followed by expression of
#' genes in that pathway in heatmap format.
#'
#' @param cutoff If the `gene_list` argument is longer than a cutoff length,
#' then the `gene_list` will be truncated into lists of genes, each of which is
#' no longer than the cutoff length. Then each heatmap will be generated for
#' each sublist
pathway_heat_pca <- function (x, gene_list, save_dir, cutoff=70, 
                              highlight_cells=NULL, group.by=c('revised', 'date'),
                              DE_dir=NULL, DE_label=NULL){
        x <- Seurat::ScaleData (x, features=gene_list)
        y <- Seurat::RunPCA (x, features=gene_list)
        pdf (save_dir, width=20, height=10)
        if (is.null (highlight_cells)){
                print (plot_dim_red (y, by_group= group.by, DR='pca', all_labels=T))
        }else{
                print (gg_plot_dim_red (y, by_group = group.by, DR='pca' ,
                                    size_highlight=highlight_cells, highlight_font=2, dims=c(1,2)))
        }
        if (!is.null (DE_dir) & !is.null (DE_label) ){
                print (plot_gene_PC (y, directory=DE_dir, label=DE_label, color_by='revised'))
        }else{print ( plot_gene_PC (y) )}
        rm (y)

        heat_high_not (x, gene_list, highlight_cells)
        if (length (gene_list) > cutoff){
                print ('making separate heatmaps')
                for (i in seq (1, length(gene_list), cutoff ) ){
                        if (i < length (gene_list) ){
                                short_list <- gene_list [i:(i+cutoff)]
                                heat_high_not (x, short_list, highlight_cells)
                        }
                }
        }
        dev.off ()
}


#' Draw the heatmap and PCA on genes of different pathways
#'
#' @param x a Seurat object
#' @param save_dir where outputs will be saved. For each pathway, one pdf file
#' will be made, including the PCA in the first page, heatmap the second, and
#' truncated heatmaps if the heatmap is too dense
#' @param all_path a dataframe with 2 columns: 'pathway' for the common names
#' of the pathway that is used for naming the output file, and 'kegg_id' for
#' the KEGG ID for that pathway. This ID may be obtained from:
#' https://www.genome.jp/kegg/pathway.html
#' @param highlight_cells which cells to highlight in the PCA and heatmap
#' @export
all_pathways_heat_pca <- function (x, save_dir, all_path=NULL,
                                   highlight_cells=NULL, DE_dir=NULL,
                                   DE_label=NULL){
        if (is.null(all_path)){data (KeggID); all_path <- KeggID
        }
        if (!dir.exists (save_dir)){dir.create (save_dir)}
        for (one_pathway in all_path$pathway){
                print (paste ('analysing', one_pathway, 'pathway') )
                one_path <- get_path_genes (all_path, one_pathway)
                if (length (one_path) !=0 ){
                        pathway_heat_pca ( x, one_path,  paste (save_dir, 
                                paste (one_pathway, '.pdf', sep=''), sep='/'),
                                highlight_cells=highlight_cells, DE_dir=DE_dir,
                                DE_label=DE_label)
                }
        }
}

#' Draw the heatmap for the eigengene for single cell
#'
#' @description Due to the problem with the sign of eigenvector, this function
#' is not recommendd to visualise pathway activity. Instead, you may use the
#' `get_module_score` function and plot the results using `seurat_heat`
#' @param x a Seurat object
#' @param threshold This function assigns the cells with the highest total gene
#' expressions in a particular pathway to have positive eigen values. This
#' `threshold` argument sets the range of cells to calculate the total
#' expression. For example, a threshold of 0.02 means that the total expression
#' will be calculated for the bottom 2% and top 2% of the cells
#' @examples
#' path_eigen (all_data)
#' @export
path_eigen <- function (x, all_path=NULL, threshold=0.02,
                        select_cells=NULL, return_mat=F){
        if (is.null(all_path)){data (KeggID); all_path <- KeggID
        }
        eigengene <- list()
        for (i in 1:length (all_path$pathway) ){
                print (paste ('analysing', all_path$pathway[i], 'pathway') )
                pgenes <- get_path_genes (all_path, all_path$pathway [i] )
                if (length (pgenes) != 0){
                        y <- ScaleData (x [pgenes, ], features=pgenes)
                        y <- RunPCA (y, features=pgenes)

                        eig_gene <- y@reductions[['pca']]@cell.embeddings [, 1]

                        # set up the sign of eigengene
                        thres <- quantile (eig_gene, c(threshold, 1.- threshold))
                        down_eig <- which (eig_gene < thres [1])
                        down_exp <- sum ( as.matrix (y [['RNA']]@data [, down_eig] ))

                        up_eig <- which (eig_gene > thres [2])
                        up_exp <- sum ( as.matrix (y [['RNA']]@data [, up_eig]) )

                        if (down_exp < up_exp){
                                eigengene [[i]] <- eig_gene
                        }else{
                                print ('flipping the sign of the eigengene')
                                eigengene [[i]] <- -eig_gene
                        }
                        rm (y)
                }else{eigengene [[i]] <- rep (NA, dim (x)[2] )}
        }

        eigengene <- do.call (cbind, eigengene)
        na_columns <- apply ( eigengene, 2, function (x) {sum (is.na (x)) > 0 } )
        colnames (eigengene) <- all_path$pathway
        eigengene <- eigengene [, !na_columns]
        metadata <- x@meta.data [  match ( rownames (eigengene), colnames (x) ), ]
        eigen_seurat <- Seurat::CreateSeuratObject ( t(eigengene), meta.data=metadata  )

        if (return_mat){return (eigen_seurat)
        }else{
                if (is.null (select_cells) ){
                        seurat_heat (eigen_seurat, color_row=rownames (eigen_seurat), group.by
                                         = c('revised', 'date'), slot='data', column_rotation=90,
                                         color_scale=c('blue', 'gray', 'red'))
                }else{
                        seurat_heat_highlight (eigen_seurat, select_cells, rownames (eigen_seurat), 
                                                   c('revised', 'date'), average=T, 
                                                   color_scale = c('blue', 'gray', 'red'))
                }
        }
}

#' Calculate the module score of KEGG pathway
#' 
#' @description It is extremely important to set seed for this function
#' @param x a Seurat object
#' @param save_path where to save to output dataframe. It is recommended to
#' save the results because this function tends to be time consuming.
#' @param all_path a dataframe with `pathway` column being the comon pathway
#' names, and the `kegg_id` column being the KEGG ID for the pathways
#' @param pgenes a named list of vectors of genes
#' @param append_meta if TRUE, a Seurat object will be created with the
#' expression matrix being the module scores and metadata inherited from the
#' input Seurat object. This output will be more convenient for doing violin
#' plot and heatmap
#' @return a matrix or a Seurat object whose rows are gene module names in
#' `pgenes` and whose columns are cell names
#' @export
get_module_score <- function (x, save_path, all_path=NULL, pgenes=NULL, append_meta=F){
        if (is.null(all_path)){data (KeggID); all_path <- KeggID
        }
        if ( file.exists (save_path) ){
                print ('loading saved module scores')
                module_scores <- utils::read.csv (save_path, row.names=1)
                # e.g. PI3K-Akt signaling may be changed to PI3K.Akt, need to
                # change it back
                rownames (module_scores) <- gsub ('\\.', '-', rownames (module_scores) )
                colnames (module_scores) <- gsub ('^X', '', colnames (module_scores))
        }else{
                if (is.null (pgenes)){
                        pgenes <- list ()
                        for (i in 1:length (all_path$pathway) ){
                                print (paste ('geting genes from', all_path$pathway[i], 'pathway') )
                                pgenes [[i]] <- get_path_genes (all_path, all_path$pathway [i] )
                        }

                        names (pgenes) <- all_path$pathway
                        # remove modules with 0 genes
                        module_num <- sapply (pgenes, length) 
                        pgenes <- pgenes [ module_num != 0 ]
                }
                names (pgenes) <- paste (names (pgenes), '_', sep='')
                # calculate module scores
                x <- Seurat::AddModuleScore(x, features = pgenes, name = names (pgenes) )
                # this function tends to append numbers after gene module names, need
                # to clean them using regexp
                module_names <- paste (names (pgenes), 1:length(pgenes), sep='')
                module_names <- gsub ('-', '.', module_names)
                module_scores <- t(x@meta.data [, module_names])
                rownames (module_scores) <- gsub ('_[0-9]+$', '', rownames (module_scores))
                utils::write.csv (module_scores, save_path)
        }

        if (!append_meta){return (module_scores)
        }else{
                metadata <- x@meta.data 
                module_scores <- module_scores [, match (rownames (metadata), colnames (module_scores) ) ]
                path_seurat <- Seurat::CreateSeuratObject (module_scores, meta.data=metadata)
                return (path_seurat)
        }
}
