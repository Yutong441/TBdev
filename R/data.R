#' Configuration of aesthetic parameters
#'
#' @format a list
#' \describe{
#'         \item{color_vec}{a named vector of colors to plot for specific names}
#'         \item{remove_keys}{what terms to remove from GO/KEGG/Reactome}
#'         \item{fontsize}{the font size of all graphical elements except
#'         `ggplot2::geom_text` and `ggrepel::geom_text_repel`}
#'         \item{point_fontsize}{the font size of labelled text in the graph by
#'         `geom_text` and `geom_text_repel`}
#'         \item{legend_point_size}{the point size in legend symbol, if the value is
#'         discrete}
#'         \item{normal_shape}{shape of points, must be between 21 and 24}
#'         \item{highlight_shape}{shape of highlighted points in `plot_dim_red
#'         (size_highlight=)` argument, must be between 21 and 24}
#'         \item{font_fam}{font family}
#'         \item{gfont_fam}{for `grid::gpar` argument, 'sans' mans 'Arial'}
#'         \item{highlight_font}{a list containing the fontface and fontsize. This is
#'         for the labelling letters in arranging multiple plots using `arrange_plots`}
#'         \item{ridge_alpha}{transparency of ridge plot}
#'         \item{heatmap_color}{a vector of 3 items corresponding to the color of low,
#'         middle and high range of the heatmap}
#'         \item{date_color_vec}{a vector for plotting points related to the `date`
#'         metadata field in a Seurat object}
#'         \item{palette}{color palette to use for `ggplot2::scale_color_continuous`}
#'         \item{arrow_angle}{the angle of arrow when plotting the arrow axis}
#'         \item{arrow_length}{the length of the head of the arrow}
#'         \item{arrow_length_unit}{the measurement unit of `arrow_length`}
#'         \item{arrow_type}{the shape of the arrow head}
#'         \item{arrow_thickness}{the thickness of the arrow}
#'         \item{arrow_linejoin}{the shape of the lateral ends of the arrow head}
#'         \item{cell_order}{the order of items that appear in the legend}
#' }
"format_conf"

#' Simplify the GO/KEGG/Reactome terms
#'
#' @format a dataframe
#' \describe{
#'         \item{ori}{the term to be replaced}
#'         \item{sub}{the simplified term}
#' }
"GOsimp"

#' A list of transcriptional factors
#'
#' @format a vector
#' @references
#' \url{https://www.cell.com/cell/comments/S0092-8674(18)30106-5#supplementaryMaterial}
#' \url{http://humantfs.ccbr.utoronto.ca/download.php}
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
#' \describe{
#'         \item{pathway}{the common name of a pathway}
#'         \item{kegg_id}{the KEGG ID of that pathway}
#' }
"KeggID"

#' Cell type collection
#'
#' @format a list of cell types
#' \describe{
#'         \item{TB_lineage}{the cell types within the trophoblast lineage}
#'         \item{non_emb_lineage}{the cell types not within the embryonic lineage}
#'         \item{in_vitro_cells}{the in vitro cell types}
#' }
"CT"

#' Immunofluorescence quantification data
#'
#' @format a dataframe with the following columns for each cell (along the
#' rows) at its maximum intensity projection:
#' \describe{
#'        \item{condition}{the signaling pathway tested}
#'        \item{major_axis_length.struct}{The major axis length of an entire
#'        cell, assuming cell is an ellipse}
#'        \item{minor_axis_length.struct}{The minor axis length of an entire
#'        cell, assuming cell is an ellipse}
#'        \item{perimeter.struct}{the perimeter of an entire cell}
#'        \item{area.struct}{the area of an entire cell}
#'        \item{major_axis_length}{The major axis length of a nucleus
#'        nucleus, assuming a nucleus is an ellipse}
#'        \item{minor_axis_length}{The major minor axis length of a nucleus
#'        nucleus, assuming a nucleus is an ellipse}
#'        \item{perimeter}{the perimeter of a nucleus}
#'        \item{area}{the area of a nucleus}
#'        \item{series}{The series number in the original .lif file used for
#'        blinding}
#'        \item{HLAG}{The fluorescence level of HLA-G in a cell}
#'        \item{CGB}{The fluorescence level of CGB in a cell}
#'        \item{TFAP2C}{The fluorescence level of TFAP2C in the nucleus of a cell}
#'        \item{HLAG_positive}{whether the cell is HLA-G positive}
#'        \item{CGB_positive}{whether the cell is CGB positive}
#'        \item{TFAP2C_positive}{whether the cell is TFAP2C positive}
#'        \item{nuc_circ}{The nuclear circularity}
#'}
#' NB: 
#' 1. due to the large size of the file, it is not included in the github
#' package. Please download it following the README.md file instruction.
#' 2. All units of measurements of cell morphology are in pixel number.
#' 3. All units of measurements of fluorescence are scaled from 0 to 255 per
#' volume.
#' "all_cleaned"
