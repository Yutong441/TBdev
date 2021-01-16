# I have to copy a lot of source codes here because I wish to modify the names
# shown in the piechart. The original code does not allow me to do this.
# Because the object is ggraph instead of a simple dataframe, I cannot do just
# simply add another geom_text layer controlling the text labels.
# Reference: https://rdrr.io/bioc/enrichplot/api/

update_n <- function(x, showCategory) {
    if (!is.numeric(showCategory)) {
        return(showCategory)
    }

    ## geneSets <- geneInCategory(x) ## use core gene for gsea result
    n <- showCategory
    if (nrow(x) < n) {
        n <- nrow(x)
    }

    return(n)
}
merge_compareClusterResult <- function(yy) {
    yy_union<- yy[!duplicated(yy$ID),]
    yy_ids <- lapply(split(yy, yy$ID), function(x) {
        ids <- unique(unlist(strsplit(x$geneID, "/")))
        cnt <- length(ids)
        list(ID=paste0(ids, collapse="/"), cnt=cnt)
    })

    ids <- vapply(yy_ids, function(x) x$ID, character(1))
    cnt <- vapply(yy_ids, function(x) x$cnt, numeric(1))

    yy_union$geneID <- ids[yy_union$ID]
    yy_union$Count <- cnt[yy_union$ID]
    yy_union$Cluster <- NULL
    yy_union
}

get_y_union <- function(y, showCategory){
    y_union <- merge_compareClusterResult(y)

    n <- update_n(y_union, showCategory)
    if (is.numeric(n)) {
        y_union <- y_union[1:n,]
    } else {
        y_union <- y_union[match(n, y_union$Description),]
        n <- length(n)
    }
     if (n == 0) {
        stop("no enriched term found...")
    }

   return(y_union)

}

get_igraph <- function(x, y,  n, color, cex_line, min_edge){
    geneSets <- DOSE::geneInCategory(x) ## use core gene for gsea result
    if (is.numeric(n)) {
        y <- y[1:n, ]
    } else {
        y <- y[match(n, y$Description),]
        n <- length(n)
    }

    if (n == 0) {
        stop("no enriched term found...")
    }
    
    g <- emap_graph_build(y = y, geneSets = geneSets, color = color,
             cex_line = cex_line, min_edge = min_edge,
             pair_sim = x@termsim, method = x@method)
}

#' @importFrom igraph E "E<-"
#' @importFrom igraph V "V<-"
#' @noRd
emap_graph_build <- function(y, geneSets, color, cex_line, min_edge, 
                             pair_sim  = NULL, method = NULL) {

    if (!is.numeric(min_edge) | min_edge < 0 | min_edge > 1) {
    	stop('"min_edge" should be a number between 0 and 1.')
    }

    if (is.null(dim(y)) | nrow(y) == 1) {  # when just one node
        g <- igraph::graph.empty(0, directed=FALSE)
        g <- igraph::add_vertices(g, nv = 1)
        V(g)$name <- as.character(y$Description)
        V(g)$color <- "red"
    } else {
        w <- pair_sim
        if (method == "JC") {
            w <- w[as.character(y$Description), as.character(y$Description)]
        } else {
            w <- w[y$ID, y$ID]
        }
    }

    wd <- reshape2::melt(w)
    wd <- wd[wd[,1] != wd[,2],]
    # remove NA
    wd <- wd[!is.na(wd[,3]),]
    if (method != "JC") {
        # map id to names
        wd[, 1] <- y[wd[, 1], "Description"]
        wd[, 2] <- y[wd[, 2], "Description"]
    }

    g <- igraph::graph.data.frame(wd[, -3], directed=FALSE)
    E(g)$width <- sqrt(wd[, 3] * 5) * cex_line

    # Use similarity as the weight(length) of an edge
    E(g)$weight <- wd[, 3]
    g <- igraph::delete.edges(g, E(g)[wd[, 3] < min_edge])
    idx <- unlist(sapply(V(g)$name, function(x) which(x == y$Description)))
    cnt <- sapply(geneSets[idx], length)
    V(g)$size <- cnt
    colVar <- y[idx, color]
    V(g)$color <- colVar
    return(g)
}

prepare_pie_category <- function(y, pie = "equal") {
    pie <- match.arg(pie, c("equal", "count", "Count"))
    if (pie == "count") pie <- "Count"

    pie_data <- y[,c("Cluster", "Description", "Count")]
    pie_data[,"Description"] <- as.character(pie_data[,"Description"])
    prepare_pie_data(pie_data, pie = pie)
}

prepare_pie_data <- function(pie_data, pie = "equal",type = "category") {
    if(type == "category"){
        ID_unique <- unique(pie_data[,2])
    } else {
        ID_unique <- unique(pie_data[,3])
    }

    Cluster_unique <- unique(pie_data[,1])
    ID_Cluster_mat <- matrix(0, nrow = length(ID_unique), ncol = length(Cluster_unique))
    rownames(ID_Cluster_mat) <- ID_unique
    colnames(ID_Cluster_mat) <- Cluster_unique
    ID_Cluster_mat <- as.data.frame(ID_Cluster_mat, stringAsFactors = FALSE)
    if(pie == "Count") {
        for(i in seq_len(nrow(pie_data))) {
            ID_Cluster_mat[pie_data[i,2],pie_data[i,1]] <- pie_data[i,3]
        }
        for(kk in seq_len(ncol(ID_Cluster_mat))) {
            ID_Cluster_mat[,kk] <- as.numeric(ID_Cluster_mat[,kk])
        }
        return(ID_Cluster_mat)
    }
    for(i in seq_len(nrow(pie_data))) {
        if(type == "category"){
            ID_Cluster_mat[pie_data[i,2],pie_data[i,1]] <- 1
        } else {
            ID_Cluster_mat[pie_data[i,3],pie_data[i,1]] <- 1
    }

    }
    return(ID_Cluster_mat)
}

#' @importFrom ggplot2 aes_
#' @noRd
get_p <- function(y, g, y_union, cex_category, pie, layout){
    ## when y just have one line
    if(is.null(dim(y)) | nrow(y) == 1) {
        title <- y$Cluster
        p <- ggraph::ggraph(g) + ggraph::geom_node_point(color="red", size=5 * cex_category) +
            ggraph::geom_node_text(aes_(label=~name)) + ggplot2::theme_void() +
            ggplot2::labs(title=title)
        return(p)
    }

    if(is.null(dim(y_union)) | nrow(y_union) == 1) {
        ##return(ggraph(g) + geom_node_point(color="red", size=5) + geom_node_text(aes_(label=~name)))
        p <- ggraph::ggraph(g)
        ID_Cluster_mat <- prepare_pie_category(y, pie=pie)

        ID_Cluster_mat <- cbind(ID_Cluster_mat,1,1,0.1*cex_category)
        colnames(ID_Cluster_mat) <- c(colnames(ID_Cluster_mat)[1:(ncol(ID_Cluster_mat)-3)],
            "x", "y", "radius")


        p <- p + scatterpie::geom_scatterpie(aes_(x=~x,y=~y,r=~radius), data=ID_Cluster_mat,
                cols=names(ID_Cluster_mat)[1:(ncol(ID_Cluster_mat)-3)],
                color=NA)+
            ggplot2::xlim(-3,3) + ggplot2::ylim(-3,3) + ggplot2::coord_equal()+
            ggraph::geom_node_text(aes_(label=~name), repel=TRUE) +
            ggplot2::theme_void()+ ggplot2::labs(fill = "Cluster")
        return(p)

    }
    ggraph::ggraph(g, layout=layout)
}

#' Plot enrichment pie chart with simplified terms
#'
#' @param split separate result by 'category' variable
#' @param pie proportion of clusters in the pie chart, one of 'equal' (default) or 'Count'
#' @param legend_n number of circle in legend
#' @param pie_scale scale of pie chart or point, this parameter has been changed to "node_scale"
#' @param cex_line scale of line width
#' @param min_edge minimum percentage of overlap genes to display the edge,
#' should between 0 and 1, default value is 0.2
#'
#' @param vert_just vertically adjust the label position to prevent overlap
#' with the pie
#' @param repel_text whether to use ggrepel
#' @param force_repel how far the labels are repelled from the original
#' locations
#' @importFrom ggplot2 aes_
#' @importFrom igraph E "E<-"
#' @importFrom igraph V "V<-"
#' @author Guangchung Yu
sim_emap <- function(x, showCategory = 30, rename_vec=NULL, color = "p.adjust",
                     layout = "nicely", split=NULL, pie = "equal", legend_n =5,
                     cex_category = NULL, pie_scale = NULL, cex_line = 1,
                     min_edge=0.2, cex_label_category  = 1,
                     node_label_size = NULL, show.legend=F, vert_just=1,
                     repel_text=T, force_repel=1, AP=NULL
                     ) {
    AP <- return_aes_param (AP)
    if (!is.null(node_label_size))
        message("node_label_size parameter has been changed to 'cex_label_category'")
    # if (is.null(cex_label_category)) {
        # if (!is.null(node_label_size)) {
            # cex_label_category <- node_label_size
        # } else {
            # cex_label_category <- 3
        # }
    # }

    if (!is.null(pie_scale))
        message("pie_scale parameter has been changed to 'cex_category'")

    if (is.null(cex_category)) {
        if (!is.null(pie_scale)) {
            cex_category <- pie_scale
        } else {
            cex_category <- 1
        }
    }

    label_category <- 3
    ## pretreatment of x, just like dotplot do
    y <- ggplot2::fortify(x, showCategory=showCategory,
                                      includeAll=TRUE, split=split)
    y$Cluster <- sub("\n.*", "", y$Cluster)
    ## geneSets <- geneInCategory(x) ## use core gene for gsea result

    ## Data structure transformation, combining the same ID (Description) genes
    y_union <- get_y_union(y = y, showCategory = showCategory)
    y <- y[y$ID %in% y_union$ID, ]

    geneSets <- stats::setNames(strsplit(as.character(y_union$geneID), "/",
                                  fixed = TRUE), y_union$ID)

    g <- emap_graph_build(y=y_union,geneSets=geneSets,color=color,
                          cex_line=cex_line, min_edge=min_edge,
                          pair_sim = x@termsim, method = x@method)
    ID_Cluster_mat <- prepare_pie_category(y,pie=pie)

    # I only added these 3 lines
    if (!is.null (rename_vec)){
            V(g)$new_name <- rename_vec [match (V(g)$name,names (rename_vec) )]
            #y$Description <- rename_vec [match (y$Description,names (rename_vec) )]
            #ori_y_union <- y_union$Description
            #y_union$Description <- rename_vec [match (y_union$Description,names (rename_vec) )]
    }else{V(g)$new_name <- V(g)$name}

    p <- get_p(y = y, g = g, y_union = y_union, cex_category = cex_category,
               pie = pie, layout = layout)
    if (is.null(dim(y)) | nrow(y) == 1 | is.null(dim(y_union)) | nrow(y_union) == 1)
        return(p)


    p <- ggraph::ggraph(g, layout=layout)
    if (length(E(g)$width) > 0) {
        p <- p + ggraph::geom_edge_link(alpha=.8, aes_(width=~I(width)),
                                colour='darkgrey')
    }

    ## then add the pie plot
    ## Get the matrix data for the pie plot


    # plot the edge
    # get the X-coordinate and y-coordinate of pies
    aa <- p$data

    desc <- y_union$Description[match(rownames(ID_Cluster_mat),
                                      y_union$Description)]
    i <- match(desc, aa$name)

    ID_Cluster_mat$x <- aa$x[i]
    ID_Cluster_mat$y <- aa$y[i]

    #Change the radius value to fit the pie plot
    radius <- NULL
    ID_Cluster_mat$radius <- sqrt(aa$size[i] / sum(aa$size) * cex_category)
    #ID_Cluster_mat$radius <- sqrt(aa$size / pi)

    x_loc1 <- min(ID_Cluster_mat$x)
    y_loc1 <- min(ID_Cluster_mat$y)
    ## x_loc2 <- min(ID_Cluster_mat$x)
    ## y_loc2 <- min(ID_Cluster_mat$y)+0.1*(max(ID_Cluster_mat$y)-min(ID_Cluster_mat$y))
    if(ncol(ID_Cluster_mat) > 4) {
        p <- p + scatterpie::geom_scatterpie(aes_(x=~x,y=~y,r=~radius), data=ID_Cluster_mat,
            cols=colnames(ID_Cluster_mat)[1:(ncol(ID_Cluster_mat)-3)],color=NA) +
            coord_equal()
        if (utils::packageVersion("ggrepel") >= "0.9.0") {
            p <- p + ggraph::geom_node_text(aes_(label=~new_name), repel=repel_text,
                size = AP$point_fontsize, bg.color = NA, vjust=vert_just, force=force_repel)
        } else {
            p <- p + ggraph::geom_node_text(aes_(label=~new_name), repel=repel_text,
                       size = AP$point_fontsize, vjust=vert_just, force=force_repel, segment.color=NA)
        }
        p <- p + theme_void() 
        if (show.legend){
            p <- p+ scatterpie::geom_scatterpie_legend (
                ID_Cluster_mat$radius, x=x_loc1, y=y_loc1,
                n = legend_n,
                labeller=function(x) {round(sum(aa$size) * x^2 / cex_category) }
                ) + labs(fill = "Cluster")
        }
        return(p)
    }
    ## annotate("text", label = "gene number", x = x_loc2, y = y_loc2, size = 4, colour = "red")
    title <- colnames(ID_Cluster_mat)[1]
    p + geom_node_point(aes_(color=~color, size=~size))
    if (utils::packageVersion("ggrepel") >= "0.9.0") {
        p <- p + ggraph::geom_node_text(aes_(label=~new_name), repel=repel_text,
            bg.color = "white",
            size=AP$point_fontsize, vjust=vert_just, force=force_repel)
    } else {
        p <- p + ggraph::geom_node_text(aes_(label=~new_name), repel=repel_text,
            size=AP$point_fontsize, vjust=vert_just, force=force_repel, segment.color=NA)
    }
    p + ggplot2::theme_void() +
        ggplot2::scale_color_continuous(low="red", high="blue", name = color,
                               guide=ggplot2::guide_colorbar(reverse=TRUE)) +
        ggplot2::scale_size(range=c(3, 8) * cex_category)  +labs(title= title)
}
