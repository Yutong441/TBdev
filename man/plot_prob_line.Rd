% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/similarity.R
\name{plot_prob_line}
\alias{plot_prob_line}
\title{Probability of fitting in river plot}
\usage{
plot_prob_line(
  prob_data,
  select_cells,
  selected_lines = NULL,
  time_ind = "x",
  band_thick = 5,
  band_trans = 0.5,
  branch_ind = "branch",
  sel_branch = "branch1",
  meta = NULL,
  meta_time = "MGP_PT",
  meta_type = "broad_type",
  vjust = 1.2,
  thickness = 0.1,
  normalize_data = F,
  remove_underscore = T,
  AP = NULL
)
}
\arguments{
\item{select_cells}{a character vector corresponding to the names of the
cells to plot, the column storing in names should be indicated in \code{meta_type}}

\item{selected_lines}{plot the vertical lines intersecting at the maximal
probability for which cell types. If NULL, all vertical lines will be shown.
To turn off this feature, input NA}

\item{time_ind}{which column stores the pseudotime information for
\code{plot_data}}

\item{band_thick}{thickness of the lines}

\item{band_trans}{transparency of the lines}

\item{branch_ind}{which branch to plot}

\item{sel_branch}{column of \code{plot_data} that contains the branch information}

\item{meta}{dataframe for plotting a color bar that serves as a reference of
the cell types that occur in each pseudotime}

\item{meta_time}{which column in \code{meta} that stores the pseudotime
informaion}

\item{meta_type}{which column in \code{meta} that stores the cell type
information}

\item{vjust}{how far down the colorbar should be from the minimum point of
similarity}

\item{thickness}{thickness of the colorbar}

\item{normalize_data}{whether to normalize the probability within a given
cell type}

\item{plot_data}{dataframe for the river plot}
}
\description{
Probability of fitting in river plot
}
