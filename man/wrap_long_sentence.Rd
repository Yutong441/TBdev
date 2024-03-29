% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/line_break.R
\name{wrap_long_sentence}
\alias{wrap_long_sentence}
\title{Line break after very long words}
\usage{
wrap_long_sentence(vec, index_sepa = ",", break_sepa = "\\n", thres = 8)
}
\arguments{
\item{vec}{a character vector, e.g. \code{c('HLA-G','TBX3','CGA')}}

\item{index_sepa}{where may a line break be inserted}

\item{break_sepa}{the separator for line break, i.e. '\n' for most
applications except \verb{<br>} for markdown syntax.}

\item{thres}{how many letters maximum in each line}
}
\description{
I could have used \code{stringr::str_wrap} but that function does
not consider if the input string already has breaks
}
