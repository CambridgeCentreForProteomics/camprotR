#' Adds new feature describing the average reporter Signal/Noise ratio.
#'
#' @description PD column `Average.Reporter.SN` is NA when any tags have missing
#' values. This function adds an average reporter SN column which ignores missing
#' values. The assumption is that the intensity values in the `exprs` matrix
#' are Signal/Noise ratios, as reported by PD by default
#'
#' @param obj `MSnSet`. Contains PSMs for TMT quantification
#' @param sn_col `string`. Name of output column containing the average signal:noise.
#'
#' @return Returns an `MSnSet` with the average SN included as a new feature column.
#' @export
update_average_sn <- function(obj,
                              sn_col='Average.Reporter.SN'){

  fData(obj)[[sn_col]] <- rowMeans(exprs(obj), na.rm=TRUE)
  obj
}


#' Filter a PSM-level MSnSet to remove low quality PSMs
#'
#' @description Filter PSMs from TMT quantification to remove the following:
#'
#' 1. Missing values (`NA`) for all tags
#' 2. Interference/co-isolation above a set value (default=100, e.g no filtering)
#' 3. Signal:noise ratio below a set value (default=0, e.g no filtering)
#'
#' @param obj `MSnSet`. Contains PSMs.
#' @param inter_thresh `numeric`. Maximum allowed interference/co-isolation.
#' @param sn_thresh `numeric`. Minimum allowed signal:noise threshold.
#' @param master_protein_col `string`. Name of column containing master
#' proteins.
#' @param inter_col `string`. Name of column containing the interference value.
#' @param sn_col `string`. Name of column containing the signal:noise value.
#' @param verbose `boolean`. Default is TRUE, use verbose output messages.
#'
#' @return Returns an `MSnSet` with the filtered PSMs.
#' @export
filter_TMT_PSMs <- function(obj,
                            inter_thresh=100,
                            sn_thresh=0,
                            master_protein_col='Master.Protein.Accessions',
                            inter_col='Isolation.Interference.in.Percent',
                            sn_col='Average.Reporter.SN',
                            verbose=TRUE){

  if(verbose) message("Filtering PSMs...")

  obj <- obj[rowSums(is.finite(exprs(obj)))>0,]
  if(verbose) message_parse(fData(obj), master_protein_col, "No quant filtering")

  obj <- obj[fData(obj)[[inter_col]]<=inter_thresh,]
  if(verbose) message_parse(fData(obj), master_protein_col, "Co-isolation filtering")

  obj <- obj[fData(obj)[[sn_col]]>=sn_thresh,]
  if(verbose) message_parse(fData(obj), master_protein_col, "S:N ratio filtering")

  obj
}
