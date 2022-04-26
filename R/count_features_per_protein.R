#' Count the number of features (rows) per protein
#'
#' @description Given an MSnSet, determine the number of features (rows) per
#' protein for each sample (column).
#'
#' Useful for bias of variance correction when performing differential testing with
#' \url{https://www.bioconductor.org/packages/release/bioc/html/DEqMS.html}{DEqMS},
#' as proteins quantified by more PSMs/peptides tend to lower variance.
#'
#' @param obj `MSnSet`, rows are PSMs or peptides.
#' @param master_prot_col `string`, column name for master protein accessions.
#' Default is `"Master.Protein.Accessions"`.
#'
#' @return Returns a `tibble` where each row is a master protein, with three columns:
#'
#' 1. `sample` = sample name (i.e. column name)
#' 1. `master_protein_col` = master protein accession
#' 1. `n` = number of features assigned to the master protein
#'
#' @examples
#' # Use a small example TMT dataset included with the camprotR package
#' df <- psm_tmt_total
#'
#' # Make an MSnSet
#' df_exprs <- as.matrix(df[, grep("Abundance", colnames(df))])
#' colnames(df_exprs) <- gsub("Abundance\\.", "", colnames(df_exprs))
#'
#' df_fData <- df[, grep("Abundance", colnames(df), invert = TRUE)]
#'
#' psm <- MSnbase::MSnSet(exprs = df_exprs, fData = df_fData)
#'
#' # Count the number of PSMs per protein
#' count_features_per_protein(psm, master_prot_col = "Master.Protein.Accessions")
#'
#' @export
count_features_per_protein <- function(obj, master_prot_col = 'Master.Protein.Accessions') {
  obj %>%
    exprs() %>%
    data.frame() %>%
    tibble::rownames_to_column('feature') %>%
    tidyr::pivot_longer(-.data$feature, names_to = 'sample') %>%
    dplyr::mutate(sample = remove_x(.data$sample)) %>%
    filter(is.finite(.data$value)) %>%
    merge(fData(obj)[, master_prot_col, drop = FALSE], by.x = 'feature', by.y = 'row.names') %>%
    group_by(.data$sample, !!sym(master_prot_col)) %>%
    tally()
}
