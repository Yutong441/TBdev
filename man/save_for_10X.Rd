% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/expr_mat.R
\name{save_for_10X}
\alias{save_for_10X}
\title{Store 10X data}
\usage{
save_for_10X(x, save_dir, batch_size = 2000, gene_name = NULL, gene_col = NULL)
}
\arguments{
\item{x}{a data.table with the column names being the cell IDs}

\item{gene_name}{a list of genes for the expression matrix}

\item{gene_col}{a column in the expression matrix}
}
\description{
Store the expression matrix as a sparse matrix for easy loading
with the \code{read10X} function in Seurat
}
