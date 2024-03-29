% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dim_red.R
\name{run_dim_red}
\alias{run_dim_red}
\title{Run DR}
\usage{
run_dim_red(
  x,
  assay = NULL,
  select_cells = NULL,
  feature_num = 2000,
  find_var_features = FALSE,
  save_mem = NULL,
  normalize = FALSE,
  var_scale = FALSE,
  pca_features = NULL,
  run_diff_map = FALSE,
  run_umap = TRUE,
  cluster = TRUE,
  cluster_res = 0.5,
  neighbor_dim = 10
)
}
\arguments{
\item{assay}{assay to use for dimensionality reduction}

\item{select_cells}{run PCA on a selected population of cells;
should be a boolean vector}

\item{feature_num}{how many variably expressed genes to compute}

\item{find_var_features}{whether to run variable gene analysis}

\item{save_mem}{for large matrices whether the entire matrix should be
stored, or simply storing the \code{VariableFeatures}}

\item{normalize}{whether to normalize the data}

\item{var_scale}{the variable features only to reduce memory}

\item{pca_features}{run PCA on a selected list of genes, should be a
character vector}

\item{run_diff_map}{whether to run diffusion map using the package
destiny. Depending on the dataset, this may take several minutes}

\item{run_umap}{whether to run umap}

\item{cluster}{whether to perform clustering as well}
}
\value{
a Seurat object
}
\description{
Run DR
}
