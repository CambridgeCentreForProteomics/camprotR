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
silac_psm_seq_int <- function(obj,sequence_col='Annotated.Sequence'){
  
  obj <- obj %>%
    filter(Quan.Channel!='') %>%
    rowwise() %>%
    filter(is.finite(Precursor.Abundance)) %>%
    mutate(sequence_col=make_annotated_upper(!!sym(sequence_col))) %>%
    remove_silac_modifications(level='psm') %>%
    mutate(Modifications=psm_to_peptide_style_modifications(Modifications))
  
  obj_seq <- obj %>%
    group_by(!!sym(sequence_col), Quan.Channel, Modifications) %>%
    summarise(sequenced=any(is.finite(Precursor.Abundance))) %>%
    mutate(Quan.Channel=paste0('Sequenced_', Quan.Channel)) %>%
    spread(key=Quan.Channel, value=sequenced, fill=FALSE)
  
  obj_int <- obj %>%
    group_by(!!sym(sequence_col), Modifications) %>%
    summarise(max_interference=max(Isolation.Interference.in.Percent))
  
  obj_seq_int<- merge(obj_seq, obj_int, by=c(sequence_col, 'Modifications'))
  
  return(obj_seq_int)
}
