
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
silac_psm_seq_int <- function(obj, sequence_col='Sequence'){
  
  obj <- obj %>%
    filter(Quan.Channel!='') %>%
    rowwise() %>%
    filter(is.finite(Precursor.Abundance)) %>%
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


#' Determine the SILAC incorporation rate given intensity values for light and heavy
#' 
#' @description SILAC incorporation is estimated by H / (H + L), where H=Heavy
#' and L=Light labelled peptide intensity. If L is NA, incorporation==1, if H is
#' NA, incorproration==0
#'
#' @param light `numeric` Light peptide intensity
#' @param heavy `numeric` Heavy peptide intensity
#' @return `incorporation`
#' @export
get_incorporation <- function(light, heavy){
  
  if (is.na(light)|is.na(heavy)){
    
    if (is.na(light) & is.na(heavy)){
      incorporation = NA
    }
    else if (is.na(heavy)){
      incorporation = 0
    }
    else {
      incorporation = 1
    }
  }
  else{
    incorporation = heavy/(light+heavy)
  }
  return(incorporation)
}


#' Plot annotated histogram of incorporation values
#'  
#' @description Incorporation is estimated from the set of observed features. If
#' the test sample contains cells grown in 'heavy' media only, incorporation is
#' simply the mean observed incorporation. If the test sample contains a mixture of 
#' 'Heavy' and 'Light', the incorporation is best calculated from the median
#' incorporation using the corrected 'Light' intensity. 
#' From a `data.frame`, with 
#'
#' @param obj `data.frame` containing features with column 'incorporation'.
#' If mix>0, must also contain column 'corrected_incorporation' where incorporation
#' is calculated using the corrected light peptide intensity
#' @param level `string` Name for feature level. e.g 'Peptide' or 'Protein'
#' @return `list` with `p`=`ggplot` plot and `incorporation_estimates`=`list` of 
#' incorporation estimates
#' @export
plot_incorporation <- function(obj, level='Peptide', mix=0){
  
  p <- obj %>%
    ggplot() +
    geom_histogram(aes(100*incorporation), fill='grey') +
    theme_camprot(base_size=10) +
    xlab('Observed incorporation (%)') +
    ylab('Count')
  
  total_count <- nrow(obj)
  
  if(mix>0){
    
    median_incorporation <- obj %>%
      pull(incorporation) %>%
      median(na.rm=TRUE) %>%
      "*"(100)
    
    median_corrected_incorporation <- obj %>%
      pull(corrected_incorporation) %>%
      median(na.rm=TRUE) %>%
      "*"(100)
    
    incorporation_estimates <- list(
      'Median incorporation' = median_incorporation,
      'Median incorporation (corrected)' = median_corrected_incorporation)
    text <- paste0('%ss: %s\n',
                   'Median incorporation: %#.2f %%\n',
                   'Corrected median incorporation:%#.2f %%')
    
    sprintf_values <- list(text, level, total_count, 
                           median_incorporation,
                           median_corrected_incorporation)
    
  } else{
    median_incorporation <- obj %>%
      pull(incorporation) %>%
      median(na.rm=TRUE) %>%
      "*"(100)
    
    mean_incorporation <- obj %>%
      pull(incorporation) %>%
      mean(na.rm=TRUE) %>%
      "*"(100)
    
    
    incorporation_estimates <- list(
      'Median incorporation' = median_incorporation,
      'Mean incorporation' = mean_incorporation)
    
    text <- paste0('%ss: %s\n',
                   'Median incorporation: %.2f %%\n',
                   'Mean incorporation:%.2f %%')
    sprintf_values <- list(text, level, total_count, 
                           median_incorporation, mean_incorporation)
  }
  
  label <- do.call('sprintf', sprintf_values)
  p <- p + annotate(geom='text', x=0, y = Inf, label=label, vjust=1.5, hjust=0,
                    size=3, colour=get_cat_palette(1))
  
  invisible(list('p'=p, 'incorporation_estimates'=incorporation_estimates))  
}

#' Estimate the incorporation rate directly from PD output
#'  
#' @description SILAC incorporation can be estimated from peptide-level PD output
#' and summarised at peptide or protein level.
#' 
#' @details **Peptide sequencing/mass shift identification**
#' In SILAC, Peptide identity of ions can be established through MS2 fragmentation
#' ('peptide sequencing') or by mass shift (within tolerance limits) relative to
#' a sequenced peptide. Ideally, the correlation between heavy and light intensities
#' is the same regardless of whether light or heavy peptides are identified by
#' mass shift or sequencing. If e.g light peptide intensities are not well correlated 
#' when they are identified by mass shift, this may indicate the mass shift identifications
#' are erroneously picking up 'ghost peptide', which will make incorporation estimation
#' difficult
#' 
#' **Mixing Heavy and Light material for incorporation rate testing**
#' To get around the issue of 'ghost peptide', one can spike in Light material (
#' at cell or protein extract-level) to the Heavy material being analysed.
#  
#' 
#' @param psm_infile `string` filepath to psm-level PD output
#' @param peptide_infile `string` filepath to peptide-level PD output
#' @param crap_fasta `string`. filepath to cRAP fasta
#' @param master_protein_col `string`. Name of column containing master
#' proteins.
#' @param protein_col `string`. Name of column containing all protein
#' matches.
#' @param sequence_col `string`. Name of column containing peptide sequences
#' @param modifications_col `string`. Name of column containing peptide modifications
#' @param mix `numeric`. If Light material has been spiked in, what is the abundance
#' relative to the heavy material. E.g if they are equal, mix=1.
#' Default value = 0, e.g no Light spike in
#' @param outdir `string` filepath to directory for plots and incorporation summary table
#' @return `incorporation`
#' @export
estimate_incorporation <- function(psm_infile,
                                   peptide_infile,
                                   crap_fasta,
                                   master_protein_col="Master.Protein.Accessions",
                                   protein_col="Protein.Accessions",
                                   sequence_col='Sequence',
                                   modifications_col='Modifications',
                                   mix=0,
                                   outdir){
  
  if(!dir.exists(outdir)) dir.create(outdir)
  
  psm_sequenced_data <- silac_psm_seq_int(read.delim(psm_infile),
                                          sequence_col=sequence_col)
  
  # Load the cRAP FASTA used for the PD search
  crap.fasta <- Biostrings::fasta.index(crap_fasta, seqtype = "AA")
  
  # Extract the non cRAP UniProt accessions associated with each cRAP protein
  crap.accessions <- crap.fasta %>% 
    pull(desc) %>% 
    stringr::str_extract_all("(?<=\\|).*?(?=\\|)") %>% 
    unlist()
  
  peptide_data <- parse_features(read.delim(peptide_infile), 
                                 master_protein_col=master_protein_col,
                                 protein_col=protein_col,
                                 unique_master=FALSE,
                                 silac=TRUE,
                                 level="peptide",
                                 filter_crap=TRUE,
                                 crap_proteins=crap.accessions,
                                 filter_associated_crap=TRUE) %>%
    filter(grepl('K|R|k|r', !!sym(sequence_col))) %>%
    filter(Quan.Info!='Redundant')
  
  light_col <- colnames(peptide_data)[
    grepl('Abundances.Grouped.(F\\d*.)?Light', colnames(peptide_data))]
  heavy_col <- colnames(peptide_data)[
    grepl('Abundances.Grouped.(F\\d*.)?Heavy', colnames(peptide_data))]
  
  peptide_data <- peptide_data %>%
    select(!!sym(sequence_col), !!sym(modifications_col),
           !!sym(protein_col), !!sym(master_protein_col), Number.of.Missed.Cleavages,
           'light'=all_of(light_col), 'heavy'=all_of(heavy_col)) %>%
    mutate(corrected_light=light-(mix*heavy)) %>%
    rowwise() %>%
    mutate(incorporation=get_incorporation(light, heavy),
           corrected_incorporation=get_incorporation(corrected_light, heavy))
  
  peptide_data <- peptide_data %>%
    remove_silac_modifications(level='peptide') %>%
    merge(psm_sequenced_data, by=c(sequence_col, modifications_col)) %>%
    mutate(found_light=ifelse(Sequenced_Light, 'Light: Spectrum matched', ifelse(
      is.finite(light), 'Light: Detected by mass shift', 'Light: Not found')),
      found_heavy=ifelse(Sequenced_Heavy, 'Heavy: Spectrum matched', ifelse(
        is.finite(heavy), 'Heavy: Detected by mass shift', 'Heavy: Not found')))
  
  p_cor <- peptide_data %>%
    filter(is.finite(light), is.finite(heavy)) %>%
    ggplot(aes(log10(heavy), log10(light))) +
    geom_point(size=0.2) + 
    theme_bw() +
    theme(aspect.ratio=1) +
    facet_wrap(found_heavy~found_light) +
    xlab('Heavy intensity (log10)') +
    ylab('Light intensity (log10)') +
    xlim(3, NA) +
    ylim(3, NA)
  
  ggsave(file.path(outdir, 'heavy_light_correlation.png'), p_cor)
  
  if(mix>0){
    p_cor2 <- p_cor + geom_abline(slope=1, linetype=2,
                                  intercept=log10((mix/2)/(1-(mix/2))),
                                  colour=get_cat_palette(1))
    
    ggsave(file.path(outdir, 'heavy_light_correlation2.png'), p_cor2)
    
    peptide_incorporation <- peptide_data %>%
      filter(is.finite(light), is.finite(heavy)) %>%
      plot_incorporation(mix=mix)
    
  } else {
    peptide_incorporation <- plot_incorporation(peptide_data)
  }
  
  
  ggsave(file.path(outdir, 'peptide_incorporation.png'), peptide_incorporation$p)
  
  multi_peptide <- peptide_data %>%
    # first run unique on the protein, sequence columns so each peptide sequence is
    # only counted once
    select(!!sym(master_protein_col), Sequence) %>% # 
    unique() %>%
    group_by(!!sym(master_protein_col)) %>%
    tally() %>%
    filter(n>1)
  
  # Use MsnSets here and make other combine methods available?
  if(mix>0){
    protein_data <- peptide_data %>%
      filter(is.finite(light), is.finite(heavy)) %>%
      group_by(!!sym(master_protein_col)) %>%
      summarise(corrected_incorporation=median(corrected_incorporation, na.rm=TRUE),
                incorporation=median(incorporation, na.rm=TRUE)) %>%
      merge(multi_peptide, by=master_protein_col)
  } else{
    protein_data <- peptide_data %>%
      group_by(!!sym(master_protein_col)) %>%
      summarise(corrected_incorporation=median(corrected_incorporation, na.rm=TRUE),
                incorporation=median(incorporation, na.rm=TRUE)) %>%
      merge(multi_peptide, by=master_protein_col)
  }
  
  
  protein_incorporation <- plot_incorporation(protein_data, level='Protein', mix=mix)
  ggsave(file.path(outdir, 'protein_incorporation.png'), protein_incorporation$p)
  
  incorporation_estimates <- data.frame(matrix(
    unlist(c(peptide_incorporation$incorporation_estimates,
             protein_incorporation$incorporation_estimates)),
    ncol=length(peptide_incorporation$incorporation_estimates), byrow=TRUE)) %>%
    setNames(names(peptide_incorporation$incorporation_estimates)) %>%
    mutate('level'=c('peptide', 'protein'))#
  
  write.table(incorporation_estimates, file.path(outdir, 'incorporation.tsv'), sep='\t', row.name=FALSE)
  
  invisible(NULL)
}