% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/combine_runs.R
\name{compare_IDs}
\alias{compare_IDs}
\title{Identify the proportion of exact matches in protein IDs}
\usage{
compare_IDs(seq2new_master_prot, master_prot_col)
}
\arguments{
\item{seq2new_master_prot}{\code{data}. Output from merging
\code{\link{read_infiles}} & \code{\link{get_sequence2protein}}}

\item{master_prot_col}{\code{string} Column name for master protein(s)}
}
\value{
Prints the proportion of exact matches in the console.
}
\description{
With object containing original and updated master protein
assignments, identify the proportion of exact matches in protein IDs
}
\keyword{internal}
