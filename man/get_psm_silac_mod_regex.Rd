% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/modifications.R
\name{get_psm_silac_mod_regex}
\alias{get_psm_silac_mod_regex}
\title{Get a pre-defined regex for a SILAC modification in PSM format}
\usage{
get_psm_silac_mod_regex(silac_mod)
}
\arguments{
\item{silac_mod}{SILAC modification name (call this function without
arguments to see available values).}
}
\value{
Returns a \verb{character vector}, regex for SILAC modification.
}
\description{
This function returns a regex which can be used to remove SILAC
modifications from the modification column in the PSM-level PD output using
\link[camprotR]{remove_silac_modifications}. Call without any arguments to see
a description of the available modifications.
}
