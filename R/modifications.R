#' Transformes PSM-level Modifications annotation style to peptide-level style
#'
#' @description Proteome Discoverer has a different style of annotation for
#' peptide modifications in PSM and peptide-level outputs. In PSM-level output,
#' each modification is listed separately, with the modification in parentheses
#' after the position. In peptide-level output, the modifications are summarised
#' more succintly, with the number of each modification followed by the positions
#' For example:
#' PSM: N-Term(Prot)(Acetyl); C2(Carbamidomethyl); C16(Carbamidomethyl)
#' Peptide:1xAcetyl \[N-Term\]; 2xCarbamidomethyl \[C2; C16\]
#'
#' @param psm_style_modifications `vector` of PSM-style modification annotations
#' @return `vector` of Peptide-style modification annotations
#' @examples
#' psm_mod <- "N-Term(Prot)(Acetyl); C2(Carbamidomethyl); C16(Carbamidomethyl)"
#' target <- "1xAcetyl [N-Term]; 2xCarbamidomethyl [C2; C16]"
#' psm_to_peptide_style_modifications(psm_mod) == target
#' @export
psm_to_peptide_style_modifications <- function(psm_style_modifications) {
  elements <- strsplit(psm_style_modifications, split = '; ')

  elements %>%
    lapply(function(x) {
      mods <- gsub('\\S+\\(|)$', '', x)
      amino_acids <- gsub('\\(\\S+$', '', x)

      count_mods <- table(mods)
      output_string <- vector("list", length = length(count_mods))
      for (ix in seq_along(names(count_mods))) {
        mod = names(count_mods)[[ix]]
        mod_count <- count_mods[[mod]]
        mod_positions <- which(mods %in% mod)
        output_string[[ix]] <- sprintf('%sx%s [%s]', mod_count, mod, paste(amino_acids[mod_positions], collapse = '; '))
      }

      return(paste(output_string, collapse = '; '))
    }) %>%
    unlist()
}


#' Remove SILAC labels from modifications
#'
#' @description This function removes SILAC heavy/light labels from a PD
#' Modifications column. When dealing with PSM level Proteome Discoverer output,
#' this is required to match PSMs for the same peptide, or to transfer PSM-level
#' features to peptide-level data.
#'
#' @param mod_col `character vector` Modification column from Proteome Discoverer
#' @param level `character` Either 'psm' or 'peptide'
#' @param psm_modfication_regexes `character vector` One or more regexes to match the expected SILAC modifications
#' @return `character vector` updated Modifications column
#' @export
remove_silac_modifications <- function(mod_col, level = 'psm',
                                       psm_modfication_regexes=c(get_psm_silac_mod_regex('R_13C6_15N4'),
                                                                 get_psm_silac_mod_regex('K_13C6_15N2'))){
  if (!level %in% c('psm', 'peptide')) stop('level must be psm or peptide')

  if (level == 'psm') {

    psm_modfication_regex <- sprintf(
      '(; )?(%s)',
      paste(psm_modfication_regexes, collapse='|')
      )

    mod_col <- gsub(
      '^(; )', '',
      gsub(
        psm_modfication_regex,
        '', mod_col
      )
    )
  } else {
    mod_col <- gsub(
      '^(; )', '',
      gsub(
        '(; )?\\dxLabel:\\S+ \\[.*]',
        '', mod_col
      )
    )
  }

  return(mod_col)
}

#' Get a pre-defined regex for a SILAC modification in PSM format
#'
#' @description This function returns a regex which can be used to remove SILAC
#' modifications from the modification column in the PSM-level PD output using
#' \link[camprotR]{remove_silac_modifications}. Call without any arguments to see
#' a description of the available modifications.
#'
#' @param `silac_mod` SILAC modification name (call without arguments to see available values)
#' @return `character vector` regex for SILAC modification
#' @export
get_psm_silac_mod_regex <- function(silac_mod){

  regexes <- list('R_13C6_15N4' =
                    list('desc'='Heavy R (10)',
                         'regex'='R\\d{1,2}\\(Label:13C\\(6\\)15N\\(4\\)\\)'),
                  'K_13C6_15N2' =
                     list('desc'='Heavy K (8)',
                          'regex'='K\\d{1,2}\\(Label:13C\\(6\\)15N\\(2\\)\\)'),
                  'R_13C6' =
                    list('desc'='Median R (6)',
                         'regex'='R\\d{1,2}\\(Label:13C\\(6\\)\\)'),
                  'K_2H4' =
                    list('desc'='Median K (4)',
                         'regex'='K\\d{1,2}\\(Label:2H\\(4\\)\\)'))

  if(missing(silac_mod)){

    table <- do.call('rbind', regexes) %>%
      data.frame() %>%
      tibble::rownames_to_column('name')

    return(table)
  }
  else{
    silac_mod = match.arg(silac_mod, choices=names(regexes))
    return(regexes[[silac_mod]][['regex']])
  }
}


