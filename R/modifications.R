#' Transformes PSM-level Modifications annotation style to peptide-level style
#'
#' @description Proteome Discoverer has a different style of annotation for
#' peptide modifications in PSM and peptide-level outputs. In PSM-level output,
#' each modification is listed separately, with the modification in parentheses
#' after the position. In peptide-level output, the modifications are summarised
#' more succintly, with the number of each modification followed by the positions
#' For example:
#' PSM: N-Term(Prot)(Acetyl); C2(Carbamidomethyl); C16(Carbamidomethyl)
#' Peptide:1xAcetyl [N-Term]; 2xCarbamidomethyl [C2; C16]
#'
#' @param obj `vector` of PSM-style modification annotations
#' @return `vector` of Peptide-style modification annotations
#' @example
#' psm_mod <- "N-Term(Prot)(Acetyl); C2(Carbamidomethyl); C16(Carbamidomethyl)"
#' target <- "1xAcetyl [N-Term]; 2xCarbamidomethyl [C2; C16]"
#' psm_to_peptide_style_modifications(psm_mod)==target
#' @export
psm_to_peptide_style_modifications <- function(psm_style_modifications){
  elements <- strsplit(psm_style_modifications, split='; ')

  elements %>% lapply(function(x){
    mods <- gsub('\\S+\\(|)$', '', x)
    amino_acids <- gsub('\\(\\S+$', '', x)

    count_mods <- table(mods)
    output_string <- vector("list", length = length(count_mods))
    for(ix in seq_along(names(count_mods))){
      mod=names(count_mods)[[ix]]
      mod_count <- count_mods[[mod]]
      mod_positions <- which(mods %in% mod)
      output_string[[ix]] <- sprintf('%sx%s [%s]', mod_count, mod, paste(amino_acids[mod_positions], collapse='; '))
    }

    return(paste(output_string, collapse='; '))
  }) %>% unlist()
}


#' Remove SILAC labels from modifications
#'
#' @description This function removes SILAC heavy/light labels from the
#' Modifications column. When dealing with PSM level Proteome Discoverer output,
#' this is required to match PSMs for the same peptide, or to transfer PSM-level
#' features to peptide-level data.
#'
#' @param obj `data.frame` PSM or Peptide-level output from Proteome Discoverer
#' @param level `string` Eiher 'psm' or 'peptide'
#' @return `data.frame` with the updated Modifications column
#' @export
remove_silac_modifications <- function(obj, level='psm'){
  if(!level %in% c('psm', 'peptide')) stop('level must be psm or peptide')

  if(level=='psm'){
    obj <- obj %>%
      mutate(Modifications=gsub(
        '^(; )', '',
        gsub('(; )?(K|R)\\d{1,2}\\(Label:13C\\(6\\)15N\\((2|4)\\)\\)', '',
             Modifications)))
  } else{
    obj <- obj %>%
      mutate(Modifications=gsub(
        '^(; )', '',
        gsub('(; )?\\dxLabel:\\S+ \\[.*]', '',
             Modifications)))
  }

  return(obj)
}
