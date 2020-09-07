#' Remove features which are assigned to a protein with too few supporting
#' features in total
#'
#' @description For summarisation of PSM or peptide to protein, we need a
#' minimum number of finite values per protein per sample. Where there are two
#' few finite values for a given protein in a given sample, this function will
#' replace all values for the protein features with NA. Proteins with two few
#' finite values in all samples will be removed entirely.
#'
#' This function is useful for Label-Free Quantification but also when using
#' robust summarisation for isobaric tagging without PSM-level imputation
#'
#' @param obj `MSnSet` with PSM or peptide-level quantification
#' @param min_features `numeric` Threshold for minimum features per protein
#' @param master_protein_col `character` Column name for master protein
#' @param plot Set TRUE to plot histogram of features per protein per sample
#'
#' @return `MSnSet`
#' @export
restrict_features_per_protein <- function(
  obj,
  min_features,
  master_protein_col="Master.Protein.Accessions",
  plot=TRUE){

  feature2protein <- fData(obj) %>%
    dplyr::select(!!sym(master_protein_col)) %>%
    tibble::rownames_to_column('feature_ID')

  n_feature_per_prot <- exprs(obj) %>%
    data.frame() %>%
    tibble::rownames_to_column('feature_ID') %>%
    gather(key='sample', value='value', -feature_ID) %>%
    merge(feature2protein, by="feature_ID") %>%
    filter(is.finite(value)) %>%
    group_by(!!sym(master_protein_col), sample) %>%
    tally()

  if(plot){
    p <- ggplot(n_feature_per_prot, aes(log2(n))) +
      geom_histogram() +
      theme_camprot() +
      xlab('# features (log2)')
    print(p)
  }

  retain_mask <- fData(obj)[, master_protein_col, drop=FALSE] %>%
    tibble::rownames_to_column('feature_ID') %>%
    merge(n_feature_per_prot,
          by=master_protein_col) %>%
    mutate(retain=n>=min_features) %>%
    dplyr::select(sample, retain, feature_ID) %>%
    spread(key=sample, value=retain) %>%
    tibble::column_to_rownames('feature_ID') %>%
    as.matrix.data.frame()

  colnames(retain_mask) <- remove_x(colnames(retain_mask))

  retain_mask[is.na(retain_mask)] <- FALSE
  retain_mask <- retain_mask[rownames(obj),colnames(obj)]

  masked_exprs <- exprs(obj) * retain_mask
  masked_exprs[masked_exprs==0] <- NA

  exprs(obj) <- as.matrix(masked_exprs)

  retain_proteins <- n_feature_per_prot %>%
    filter(n>=min_features) %>%
    pull(!!sym(master_protein_col))

  return_obj <- obj[fData(obj)[[master_protein_col]] %in% retain_proteins,]

  invisible(return_obj)
}
