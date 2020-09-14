#' Extract peptide sequencing and interference information from SILAC
#' PSM-level output
#'
#' @description Proteome Discoverer does not correctly propogate from PSM to
#' peptide level output which intensities are from "sequenced" peptides, e.g MS2
#' fragmentation. This information is useful to assess the accuracy of
#' quantification in peptides identified by mass shift relative to a sequenced
#' peptide in SILAC experiments
#'
#' @param obj `data.frame` PSM-level output from Proteome Discoverer
#' @param sequence_col `string` Column with peptide sequence
#' @return `data.frame` indicating which SILAC peptides were MS2 sequenced and
#' the maximum interference across all PSMs for the peptide
#' @export
silac_psm_seq_int <- function(obj, sequence_col='Sequence'){

  obj <- obj %>%
    filter(.data$Quan.Channel!='') %>%
    rowwise() %>%
    filter(is.finite(.data$Precursor.Abundance)) %>%
    remove_silac_modifications(level='psm') %>%
    mutate(Modifications=psm_to_peptide_style_modifications(.data$Modifications))

  obj_seq <- obj %>%
    group_by(!!sym(sequence_col), .data$Quan.Channel, .data$Modifications) %>%
    summarise(sequenced=any(is.finite(.data$Precursor.Abundance))) %>%
    mutate(Quan.Channel=paste0('Sequenced_', .data$Quan.Channel)) %>%
    spread(key=.data$Quan.Channel, value=.data$sequenced, fill=FALSE)

  obj_int <- obj %>%
    group_by(!!sym(sequence_col), .data$Modifications) %>%
    summarise(max_interference=max(.data$Isolation.Interference.in.Percent))

  obj_seq_int<- merge(obj_seq, obj_int, by=c(sequence_col, 'Modifications'))

  return(obj_seq_int)
}
