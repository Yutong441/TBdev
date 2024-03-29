% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/format.R
\name{add_custom_color}
\alias{add_custom_color}
\title{Customise the colors that correspond to particular features}
\usage{
add_custom_color(feature_vec, aes_param, more_prec = 0.2, color_fill = T)
}
\arguments{
\item{feature_vec}{a vector of features for which colors will be assigned}

\item{color_fill}{whether \verb{scale_color_*} or \verb{scale_fill_*} is used}
}
\value{
a \verb{scale_*} object
}
\description{
This function works for either continuous or discrete scale
}
