# --------------------------------
# Load and merge multiple datasets
# --------------------------------

#' load all Seurat objects specified by `save_robj` in the `root` folder
#'
#' @param save_robj names of the data saved in the directory `root`
#' @param all_data a list of Seurat objects
#' NB: need to pass either `save_robj` or `all_data`
#' @export
load_all_data <- function (save_robj=NULL, root=NULL, all_data = NULL){
        if (is.null (all_data) ) {
                all_data <- list()
                num_data <- length (save_robj)
        }else{num_data <- length (all_data) }
        all_genes <- list ()

        # load data from directory
        for (i in 1:num_data){
                if (!is.null (save_robj) ){
                        r_obj <- load (paste (root, save_robj[i], sep='/'))
                        dataset <- get (r_obj)
                        rm (r_obj)
                        all_data [[i]] <- dataset
                        all_data [[i]]$paper <- strsplit (save_robj[i], '/')[[1]][[1]]
                }
                all_genes [[i]] <- rownames(all_data[[i]])
        }

        # align the gene names
        common_genes <- Reduce (intersect, all_genes)
        all_datasets <- all_data

        for (i in 1:num_data){
                all_datasets[[i]] <- all_datasets [[i]] [match (common_genes, rownames (
                                                                all_datasets[[i]])), ]
                print (dim (all_datasets[[i]]))
        }
        return (all_datasets)
}

#' Merge Seurat object
#'
#' @description Merge seurat assay data and meta.data from different Seurat
#' objects. The results will be stored in the slot 'data' and assay name 'RNA'.
#' Unlike the `merge` function provided by Seurat, this works much faster
#' because not all attributes in the Seurat objects are merged
#' @param list_obj a list of Seurat objects
#' @param assays a vector of assay names to be merged. If only one name is
#' supplied, this assay from all Seurat objects all be combined.
#' @importFrom Seurat GetAssayData
#' @export
merge_seurat <- function (list_obj, assays, slot_name='data'){
        if (length (assays) == 1){ assays <- rep (assays, length (list_obj)) }
        list_data <- list ()
        for (i in 1:length(list_obj)){
                list_obj [[i]] %>%
                        GetAssayData (slot=slot_name, assay = assays[i]) -> list_data [[i]]
                # make sure all objects have the same order of rows before `cbind`
                list_data [[i]] <- list_data [[i]] [rownames (list_data[[1]]),]
        }
        # merge assay data
        all_assays <- do.call (cbind, list_data)

        # add meta data
        all_metadata <- data.table::rbindlist (lapply (list_obj, function(x){x$ID_code <- colnames (x)
                                           return (x@meta.data) }), fill=TRUE)
        rownames (all_metadata) <- all_metadata$ID_code
        seurat_assay <- Seurat::CreateSeuratObject (all_assays, meta.data= all_metadata)
        seurat_assay[['RNA']]@data <- all_assays
        return (seurat_assay)
}

# --------------------------
# Batch correction procedure
# --------------------------

        
#' Perform Seurat batch correction
#'
#' @param data_list a list of Seurat objects
#' @export
batch_correct <- function (data_list, k_filter=NULL, dims=1:30, refer=NULL,
                           num_anchor_genes=2000){
        # NB: k.filter size should not be smaller than the smallest dataset
        all_lengths <- unlist (lapply (data_list, function(x){dim (x)[2]}))
        k_filter_max <- pmin (min(all_lengths), 200)
        if (!is.null(k_filter) ) {k_filter <- pmin (k_filter, k_flter_max)
        }else{k_filter <- k_filter_max}

        all_lengths <- unlist (lapply (data_list, function(x){dim (x)[1]}))
        num_anchor_genes <- pmin ( min(all_lengths), num_anchor_genes )
        # perform batch correction
        data_anchor <- Seurat::FindIntegrationAnchors(object.list = data_list, k.filter
                                              = k_filter, dims = dims, reference=refer,
                                                anchor.features = num_anchor_genes)
        data_integrated <- Seurat::IntegrateData(anchorset = data_anchor, dims = dims )
        Seurat::DefaultAssay (data_integrated) <- 'integrated'
        return (data_integrated)
}

#' Investigate whether varying the size of k.filter can affect the
#' result of batch correction
#'
#' @param data_list list of Seurat objects to be batch corrected
#' @param directory where the results of the batch corrections will be
#' saved
#' @param start_from the starting point of testing k.filter size. The
#' end point is 200 by default.
plot_BC_filter <- function (data_list, directory, start_from=10){
        all_lengths <- unlist (lapply (data_list, function(x){dim (x)[2]}))
        k_filter_max <- pmin (min(all_lengths), 200)
        k_filter_choice <- seq (start_from, k_filter_max, 30)

        for (k_filter in k_filter_choice){
                print (paste ('integrating datasets at filter size of', k_filter))
                data_integrated <- batch_correct (data_list, k_filter)
                data_integrated <- run_dim_red (data_integrated) # from 'dim_red.R'

                # plotting
                for (DR in c('umap', 'pca')){
                        DIR$plot_dim_red (data_integrated, DR=DR)
                        ggsave (paste (directory, paste ('batch_correct', DR,
                                                         k_filter, '.pdf', sep='_' ), sep = '/'),
                                width=12, height=7)
                }
                rm (data_integrated, data_anchor)
        }
}

# --------------
# Label transfer
# --------------

#' Transfer the known labels from selected datasets to other datasets
#'
#' @description This function uses Seurat `FindTransferAnchors` and
#' `TransferData` functions to perform label transfer
#'
#' @param merged_data a Seurat object
#' @param feature a column in the metadata of the `merged_data`, which
#' identifies the label to be transferred
#' @param by_paper select which datasets from which paper are used as a
#' reference. The default is the dataset(s) that contain missing values in the
#' field specified by the `feature` argument
#' @export
label_transfer <- function (merged_data, feature, by_paper=NULL){
        if (is.null(by_paper)){
                known_cells <- !is.na(merged_data@meta.data [, feature])
        }else{
                known_cells <- merged_data$paper %in% by_paper
        }
        query <- merged_data [, !known_cells]
        reference <- merged_data [, known_cells]
        proj_anchors <- Seurat::FindTransferAnchors(reference = reference, query = query, dims = 1:30)
        predictions <- Seurat::TransferData(anchorset = proj_anchors, refdata =
                                    reference@meta.data[,feature],  dims = 1:30)

        predicted_feature <- paste ('predicted', feature, sep='_')
        merged_data@meta.data [, predicted_feature] <- NA
        merged_data@meta.data [!known_cells, predicted_feature] <- predictions$predicted.id
        #simple copying the feature across
        merged_data@meta.data [known_cells, predicted_feature] <-
                merged_data@meta.data [known_cells, feature]
        return (merged_data)
}

plot_transfer_label <- function (merged_data, feature, DR='umap', by_paper=NULL){
        if (is.null(by_paper)){
                known_cells <- !is.na(merged_data@meta.data [, feature])
        }else{
                known_cells <- merged_data$paper %in% by_paper
        }
        selected_data <- merged_data [, !known_cells]
        predicted_feature <- paste ('predicted', feature, sep='_')
        plot1<- Seurat::DimPlot (selected_data, reduction = DR, group.by='Type', label=TRUE)
        plot2<- Seurat::DimPlot (selected_data, reduction = DR, group.by='seurat_clusters', label=TRUE)
        plot3<- Seurat::DimPlot (selected_data, reduction = DR, group.by=predicted_feature, label=TRUE)
        return (plot1 + plot2 + plot3)
}
