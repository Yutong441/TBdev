% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/GPLVM_utils.R
\name{mat_to_seurat}
\alias{mat_to_seurat}
\title{Convert matrix into a list of seurat objects based on branch assignment}
\usage{
mat_to_seurat(exp_mat, ref_meta, branch_lab, label_list)
}
\arguments{
\item{exp_mat}{a matrix with rows being the cells and columns being the
genes. There must be rownames and colnames}

\item{ref_meta}{a dataframe containing information for the cells in \code{exp_mat}}

\item{branch_lab}{a vector with the same length as the row number of
\code{exp_mat}, containing the branch assignment of each cell}

\item{label_list}{per item of returned seurat object in the return value
list, cells with which labels are chosen. It should be a list}
}
\value{
a list of seurat objects of different branch assignments
}
\description{
Convert matrix into a list of seurat objects based on branch assignment
}
