% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/expr_mat.R
\name{sparsify}
\alias{sparsify}
\title{Sparsify a large matrix}
\usage{
sparsify(x, batch_size = 2000)
}
\arguments{
\item{x}{a data.table object}

\item{batch_size}{to sparsify only a subset of the dataset at a time to
prevent memory overload. Set a lower batch_size for larger dataset.}
}
\description{
Sparsify a large matrix
}
