#' Extract peptide sequencing and interference information from SILAC
#' PSM-level output
#'
#' @description Proteome Discoverer does not correctly propagate from PSM to
#' peptide level output which intensities are from "matched" peptides, e.g MS2
#' fragmentation. This information is useful to assess the accuracy of
#' quantification in peptides identified by mass shift relative to a matched
#' peptide in SILAC experiments
#'
#' @param obj `data.frame` PSM-level output from Proteome Discoverer
#' @param sequence_col `string` Column with peptide sequence
#' @param mod_col `string` Column with modifications
#' @param include_interference `logical` Should PSM interference be included too?
#' @param interference_col `string` Column with interference/co-isolation
#' @param group_cols `string` Additional feature columns to retain, beyond
#' @param psm_modfication_regexes `character vector` One or more regexes to match the expected SILAC modifications
#' Sequence and Modification. See \link[camprotR]{get_psm_silac_mod_regex}
#' @return `data.frame` indicating which SILAC peptides were MS2 matched,
#' how many PSMs per isotope, and
#' (optionally) the maximum interference across all PSMs for the peptide
#' @export
silac_psm_seq_int <- function(
  obj,
  sequence_col='Sequence',
  mod_col='Modifications',
  include_interference=FALSE,
  interference_col='Isolation.Interference.in.Percent',
  group_cols=NULL,
  psm_modfication_regexes=c(get_psm_silac_mod_regex('R_13C6_15N4'),
                            get_psm_silac_mod_regex('K_13C6_15N2'))){

  message('camprotR::silac_psm_seq_int output has changed.
  Columns indicating whether quantification is from PSM are now prefixed with
  "matched", not "Sequenced", and tallys of PSMs per isotope are included.
  Interference is not included by default (set include_interference=TRUE)')

  # check that input files have the necessary columns
  required_cols <- c(sequence_col, mod_col,
                     "Quan.Channel", "Precursor.Abundance")

  if (!all(required_cols %in% colnames(obj))) {
    stop(
      paste("The PSM input is missing the following required columns:",
            required_cols[!required_cols %in% colnames(obj)])
    )
  }

  obj <- obj %>%
    filter(.data$Quan.Channel!='') %>%
    rowwise() %>%
    filter(is.finite(.data$Precursor.Abundance))

  obj[[mod_col]] <- remove_silac_modifications(
    obj[[mod_col]], psm_modfication_regexes=psm_modfication_regexes)

  obj[[mod_col]] <- psm_to_peptide_style_modifications(obj[[mod_col]])

  obj[[sequence_col]] <- toupper(obj[[sequence_col]])

  if(!missing(group_cols)){
    group_cols <- c(group_cols, 'Quan.Channel')
  } else{
    group_cols <- 'Quan.Channel'
  }

  obj_seq <- obj %>%
    group_by(across(all_of(c(sequence_col, group_cols, mod_col)))) %>%
    tally() %>%
    mutate('matched'=TRUE) %>%
    pivot_wider(names_from=.data$Quan.Channel,
                values_from=c(.data$n, .data$matched)) %>%
    mutate(n_Light=replace_na(.data$n_Light, 0),
           n_Heavy=replace_na(.data$n_Heavy, 0),
           matched_Light=replace_na(.data$matched_Light, FALSE),
           matched_Heavy=replace_na(.data$matched_Heavy, FALSE))

  if(include_interference){

    if (!interference_col %in% colnames(obj)) {
      stop(paste("The PSM input is missing the interference column:", interference_col))
    }

  obj_int <- obj %>%
    group_by(!!sym(sequence_col), !!sym(mod_col)) %>%
    summarise(max_interference=max(!!sym(interference_col)))

  obj_seq <- merge(obj_seq, obj_int, by=c(sequence_col, 'Modifications'))
  }


  return(obj_seq)
}

