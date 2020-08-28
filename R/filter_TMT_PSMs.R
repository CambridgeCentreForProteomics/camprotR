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
