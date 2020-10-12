#' Extract the expression matrix column medians from an MSnSet
#'
#' @param obj `MSnSet`. Contains PSMs.
#' @return `vector` of expression matrix column medians in colnames order
#' @importFrom robustbase colMedians
#' @export
get_medians <- function(obj){
  medians <- robustbase::colMedians(exprs(obj), na.rm=TRUE)
  return(medians)
}


#' Center-median normalise the expression matrix in an MSnSet using medians
#' from a reference dataset
#'
#' @description Center-median normalisation is a simple normalisation method
#' that is appropriate for relative abundance proteomics such as isobaric tagging.
#' This can be achieved with `MSnbase::normalise(method='diff.center')`. However,
#' for some experimental designs, the normalisation should be against the medians
#' in another dataset. For example, for PTM studies, one may wish to isobaric
#' tag samples, pool, and then PTM-enriched, with the enriched sample quantified
#' in a separate run to the non-enriched (total) sample. In this case, it may make
#' more sense to center-median normalise the PTM-enriched samples using the median
#' from the total samaples.
#'
#' @param obj `MSnSet`. Contains PSMs.
#' @param medians `vector, numeric`. Sample medians from reference dataset
#'
#' @return Returns an `MSnSet` with the expression matrix column center-median
#' normalised
#'
#' @export
center_normalise_to_ref <- function(obj, medians){

  medians <- medians/mean(medians)
  exprs(obj) <- t(t(exprs(obj))/medians)

  return(obj)
}
