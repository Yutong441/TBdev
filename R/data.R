#' Configuration of aesthetic parameters
#'
#' @format a list
#' \describe {
#'         \item {color_vec}{a named vector of colors to plot for specific names}
#'         \item {remove_keys}{what terms to remove from GO/KEGG/Reactome}
#'         \item {fontsize}{the font size of all graphical elements except
#'         `ggplot2::geom_text` and `ggrepel::geom_text_repel`}
#'         \item {point_fontsize}{the font size of labelled text in the graph by
#'         `geom_text` and `geom_text_repel`}
#'         \item {legend_point_size}{the point size in legend symbol, if the value is
#'         discrete}
#'         \item {normal_shape}{shape of points, must be between 21 and 24}
#'         \item {highlight_shape}{shape of highlighted points in `plot_dim_red
#'         (size_highlight=)` argument, must be between 21 and 24}
#'         \item {font_fam}{font family}
#'         \item {gfont_fam}{for `grid::gpar` argument, 'sans' mans 'Arial'}
#'         \item {highlight_font}{a list containing the fontface and fontsize. This is
#'         for the labelling letters in arranging multiple plots using `arrange_plots`}
#'         \item {ridge_alpha}{transparency of ridge plot}
#'         \item {heatmap_color}{a vector of 3 items corresponding to the color of low,
#'         middle and high range of the heatmap}
#'         \item {date_color_vec}{a vector for plotting points related to the `date`
#'         metadata field in a Seurat object}
#'         \item {palette}{color palette to use for `ggplot2::scale_color_continuous`}
#'         \item {arrow_angle}{the angle of arrow when plotting the arrow axis}
#'         \item {arrow_length}{the length of the head of the arrow}
#'         \item {arrow_length_unit}{the measurement unit of `arrow_length`}
#'         \item {arrow_type}{the shape of the arrow head}
#'         \item {arrow_thickness}{the thickness of the arrow}
#'         \item {arrow_linejoin}{the shape of the lateral ends of the arrow head}
#'         \item {cell_order}{the order of items that appear in the legend}
#' }
"format_conf"

#' Simplify the GO/KEGG/Reactome terms
#'
#' @format a dataframe
#' \describe{
#'         \item {ori}{the term to be replaced}
#'         \item {sub}{the simplified term}
#' }
"GOsimp"

#' A list of transcriptional factors
#'
#' @format a vector
#' @author Dr Penfold
"TF"

#' Genes involved in cell cycle
#'
#' @format a dataframe, containing 5 columns: G1/S, S, G2/M, M, M/G1
#' Each column contains genes characteristic of the cell cycle phase
#' @source \url{https://www.cell.com/cell/fulltext/S0092-8674%2815%2900549-8}
"cell_cycle"

#' Lineage marker genes
#' 
#' @format a named vector of genes, with the names being the lineages that they
#' mark
"lineage_markers"

#' KEGG ID conversion
#' 
#' @format a dataframe
#' \describe {
#'         \item {pathway}{the common name of a pathway}
#'         \item {kegg_id}{the KEGG ID of that pathway}
#' }
"KeggID"

#' Cell type collection
#'
#' @format a list of cell types
#' \describe{
#'         \item {TB_lineage}{the cell types within the trophoblast lineage}
#'         \item {non_emb_lineage}{the cell types not within the embryonic lineage}
#'         \item {in_vitro_cells}{the in vitro cell types}
#' }
"CT"
