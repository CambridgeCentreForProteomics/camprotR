#' Count the number of features (rows) per protein
#'
#' @description For differential testing with
#' \url{https://www.bioconductor.org/packages/release/bioc/html/DEqMS.html}{DEqMS},
#' one needs to identify the number of features, e.g PSMs or peptides per protein.
#' This function returns the number of features per protein per sample.
#'
#' @param obj `MSnSet`. Contains PSMs or peptides.
#' @param master_prot_col `character` Column name for master protein ID
#'
#' @return `tibble` with feature counts per sample per protein.
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
count_features_per_protein <- function(obj, master_prot_col='Master.Protein.Accessions'){
  obj %>%
    exprs() %>%
    data.frame() %>%
    tibble::rownames_to_column('feature') %>%
    tidyr::pivot_longer(-"feature", names_to='sample') %>%
    dplyr::mutate(sample=remove_x(.data$sample)) %>%
    filter(is.finite(.data$value)) %>%
    merge(fData(obj)[,master_prot_col,drop=FALSE], by.x='feature', by.y='row.names') %>%
    group_by(.data$sample, !!sym(master_prot_col)) %>%
    tally()
}
