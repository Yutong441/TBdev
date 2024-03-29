% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SCHeat.R
\name{seurat_heat}
\alias{seurat_heat}
\title{Heatmap for Seurat object}
\usage{
seurat_heat(
  x,
  group.by = NULL,
  color_row = NULL,
  assay = "RNA",
  slot_data = "data",
  highlight = NULL,
  highlight_names = NULL,
  show_column_names = F,
  show_row_names = T,
  show_column_anna = T,
  column_rotation = 0,
  annotation_name_side = "left",
  column_names_side = "bottom",
  column_title_side = "top",
  row_names_side = "left",
  row_title_side = "left",
  show_column_bars = T,
  row_title_fontface = "bold",
  column_title_fontface = "bold",
  cluster_columns = F,
  cluster_rows = F,
  column_split = 1,
  row_reorder_levels = NULL,
  column_reorder_levels = list(),
  reorder_column = T,
  row_scale = F,
  column_scale = F,
  center_scale = F,
  color_scale = NULL,
  default_color = F,
  quantile_val = 0.05,
  provided_color = NULL,
  AP = NULL,
  title_pos = "topleft",
  heat_name = "norm count",
  break_points = NULL,
  row_legend_labels = "DE genes",
  column_legend_labels = NULL,
  column_legend_order = NULL,
  left_HA = T,
  top_HA = T,
  row_titles = character(0),
  group_order = NULL,
  main_width = NULL,
  main_height = NULL,
  grid_height = 8,
  heat_grid_height = NULL,
  automatic = T,
  ...
)
}
\arguments{
\item{x}{Seurat object}

\item{group.by}{by which feature are the cells grouped}

\item{color_row}{a vector of gene names with the name of the vector elements
being the lineages}

\item{assay}{the Seurat assay in \code{x}}

\item{highlight}{which row names to highlght}

\item{highlight_names}{display names for the highlighted and non-highlighted
groups}

\item{show_column_names}{whether to display column names (Don't confuse it
with column title, which is group name)}

\item{show_row_names}{whether to show row names}

\item{show_column_anna}{whether to show the category for the group names for
the columns around the column color bar. There is not a similar option for
rows because this function only supports one row sidebar}

\item{column_rotation}{rotate column name labels, 90 is vertical, 0 is
horizontal, there are no other choices}

\item{annotation_name_side}{side of the annotation name for column side bar}

\item{show_column_bars}{whether to show horizontal bars for the columns. It
can be a vector of the same length as \code{group.by} to specify the group(s) to
show the horizontal bars.}

\item{cluster_columns}{whether to perform hierarchical clustering on columns}

\item{cluster_rows}{whether to perform hierarchical clustering on rows}

\item{column_split}{whether to split the columns belonging to the same group
in \code{group.by} by a thin dividing line. '1' indicate splitting. 'NA'
indicates no spliting}

\item{row_reorder_levels}{a vector indicating the order of labels for row
group names that would appear from top to bottom of the heatmap. The default
is \code{cell_order} from the list supplied to the \code{AP} argument}

\item{column_reorder_levels}{a list of vectors. Each item in the list
indicates the order of each column color bar labels. Again the default is in
the \code{AP} argument}

\item{row_scale}{whether to perform row scaling.}

\item{column_scale}{whether to perform column scaling.}

\item{center_scale}{whether to put white colors half way between the min and
max of the matrix. Otherwise, white color is placed at value 0.}

\item{color_scale}{the scale of color gradient for heatmap. By default, the
\code{heatmap_color} field of the list supplied to \code{AP} should define the color
for min, middle and max of the dataset. Otherwise, please input a function
created by \code{circlize::colorRamp2}}

\item{default_color}{whether to use the default color settings of
ComplexHeatmap for main matrix}

\item{quantile_val}{if supplied, then the color gradient would not extend
over the entire range of the matrix. The default 0.05 means that color
gradient would extend from 5\% to 95\% percentile. This is to avoid outliers
skewing the dynamic range of colors.}

\item{provided_color}{a named vector where the vector contains the colors
and the names indicate the group labels that would be colored by the
corresponding color in the sidebar.}

\item{AP}{a list of aesthetic parameters. The important fields include
\code{cell_order} for ordering group labels, \code{gfont_fam}, \code{fontsize} and
\code{heatmap_color}}

\item{title_pos}{adjustment of legend title}

\item{heat_name}{name of the main heatmap to appear in the legend}

\item{break_points}{where to show values in the heatmap gradient legend}

\item{row_legend_labels}{name of the row sidebar to appear in the legend}

\item{column_legend_labels}{name of the column sidebar in the legend. You may
supply as many items as the number of column sidebars}

\item{left_HA}{whether to show row sidebars}

\item{top_HA}{whether to show column sidebars}

\item{row_titles}{whether to show row titles. choose NULL if not}

\item{group_order}{the ordering of each column}

\item{main_width}{width of the heatmap in cm}

\item{main_height}{height of the heatmap in cm}

\item{grid_height}{legend grid height, increase this value if the spacing
between legend labels become too narrow}

\item{heat_grid_height}{grid height for the color gradient legend. If NULL,
it is the same as \code{grid_height}}

\item{automatic}{whether the default legend generation process in
ComplexHeatmap is to be used. If not, all legends will be aligned vertically
and only be arranged in separate columns if the total length exceeds the
heatmap length.}

\item{slow_data}{the slot in \code{assay}}

\item{reorder_columns}{whether to order the columns according to
\code{column_reorder_levels}}
}
\description{
Equivalent to \code{Seurat::DoHeatmap} with the addition of sidebar
along the rows (features) to group them using ComplexHeatmap. However,
currently, only one sidebar can be plotted for rows. There is no limit for
the number of side bars for columns.
}
\author{
Yutong Chen
}
