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
                                   master_protein_col = "Master.Protein.Accessions",
                                   protein_col = "Protein.Accessions",
                                   sequence_col = 'Sequence',
                                   modifications_col = 'Modifications',
                                   mix = 0,
                                   outdir) {
  if (!dir.exists(outdir)) dir.create(outdir)

  psm_sequenced_data <- silac_psm_seq_int(utils::read.delim(psm_infile),
    sequence_col = sequence_col
  )

  # Load the cRAP FASTA used for the PD search
  crap.fasta <- Biostrings::fasta.index(crap_fasta, seqtype = "AA")

  # Extract the non cRAP UniProt accessions associated with each cRAP protein
  crap.accessions <- crap.fasta %>%
    pull(desc) %>%
    regmatches(.data, gregexpr("(?<=\\|).*?(?=\\|)", .data, perl = TRUE)) %>%
    unlist()

  peptide_data <- parse_features(utils::read.delim(peptide_infile),
    master_protein_col = master_protein_col,
    protein_col = protein_col,
    unique_master = FALSE,
    silac = TRUE,
    level = "peptide",
    filter_crap = TRUE,
    crap_proteins = crap.accessions,
    filter_associated_crap = TRUE
  ) %>%
    filter(grepl('K|R|k|r', !!sym(sequence_col))) %>%
    filter(.data$Quan.Info != 'Redundant')

  light_col <- colnames(peptide_data)[
    grepl('Abundances.Grouped.(F\\d*.)?Light', colnames(peptide_data))
  ]
  heavy_col <- colnames(peptide_data)[
    grepl('Abundances.Grouped.(F\\d*.)?Heavy', colnames(peptide_data))
  ]

  peptide_data <- peptide_data %>%
    select(!!sym(sequence_col), !!sym(modifications_col),
      !!sym(protein_col), !!sym(master_protein_col), .data$Number.of.Missed.Cleavages,
      'light' = all_of(light_col), 'heavy' = all_of(heavy_col)
    ) %>%
    mutate(corrected_light = .data$light - (mix * .data$heavy)) %>%
    rowwise() %>%
    mutate(
      incorporation = get_incorporation(.data$light, .data$heavy),
      corrected_incorporation = get_incorporation(.data$corrected_light, .data$heavy)
    )

  peptide_data <- peptide_data %>%
    remove_silac_modifications(level = 'peptide') %>%
    merge(psm_sequenced_data, by = c(sequence_col, modifications_col)) %>%
    mutate(
      found_light = ifelse(.data$Sequenced_Light, 'Light: Spectrum matched', ifelse(
        is.finite(.data$light), 'Light: Detected by mass shift', 'Light: Not found'
      )),
      found_heavy = ifelse(.data$Sequenced_Heavy, 'Heavy: Spectrum matched', ifelse(
        is.finite(.data$heavy), 'Heavy: Detected by mass shift', 'Heavy: Not found'
      ))
    )

  p_cor <- peptide_data %>%
    filter(is.finite(.data$light), is.finite(.data$heavy)) %>%
    ggplot(aes(log10(.data$heavy), log10(.data$light))) +
    geom_point(size = 0.2) +
    theme_bw() +
    theme(aspect.ratio = 1) +
    facet_wrap(found_heavy ~ found_light) +
    xlab('Heavy intensity (log10)') +
    ylab('Light intensity (log10)') +
    xlim(3, NA) +
    ylim(3, NA)

  ggsave(file.path(outdir, 'heavy_light_correlation.png'), p_cor)

  if (mix > 0) {
    p_cor2 <- p_cor + geom_abline(
      slope = 1, linetype = 2,
      intercept = log10((mix / 2) / (1 - (mix / 2))),
      colour = get_cat_palette(1)
    )

    ggsave(file.path(outdir, 'heavy_light_correlation2.png'), p_cor2)

    peptide_incorporation <- peptide_data %>%
      filter(is.finite(.data$light), is.finite(.data$heavy)) %>%
      plot_incorporation(mix = mix)
  } else {
    peptide_incorporation <- plot_incorporation(peptide_data)
  }


  ggsave(file.path(outdir, 'peptide_incorporation.png'), peptide_incorporation$p)

  multi_peptide <- peptide_data %>%
    # first run unique on the protein, sequence columns so each peptide sequence is
    # only counted once
    select(!!sym(master_protein_col), .data$Sequence) %>% #
    unique() %>%
    group_by(!!sym(master_protein_col)) %>%
    tally() %>%
    filter(n > 1)

  # Use MsnSets here and make other combine methods available?
  if (mix > 0) {
    protein_data <- peptide_data %>%
      filter(is.finite(.data$light), is.finite(.data$heavy)) %>%
      group_by(!!sym(master_protein_col)) %>%
      summarise(
        corrected_incorporation = stats::median(.data$corrected_incorporation, na.rm = TRUE),
        incorporation = stats::median(.data$incorporation, na.rm = TRUE)
      ) %>%
      merge(multi_peptide, by = master_protein_col)
  } else {
    protein_data <- peptide_data %>%
      group_by(!!sym(master_protein_col)) %>%
      summarise(
        corrected_incorporation = stats::median(.data$corrected_incorporation, na.rm = TRUE),
        incorporation = stats::median(.data$incorporation, na.rm = TRUE)
      ) %>%
      merge(multi_peptide, by = master_protein_col)
  }


  protein_incorporation <- plot_incorporation(protein_data, level = 'Protein', mix = mix)
  ggsave(file.path(outdir, 'protein_incorporation.png'), protein_incorporation$p)

  incorporation_estimates <- data.frame(matrix(
    unlist(c(
      peptide_incorporation$incorporation_estimates,
      protein_incorporation$incorporation_estimates
    )),
    ncol = length(peptide_incorporation$incorporation_estimates), byrow = TRUE
  )) %>%
    setNames(names(peptide_incorporation$incorporation_estimates)) %>%
    mutate('level' = c('peptide', 'protein')) #

  utils::write.table(incorporation_estimates, file.path(outdir, 'incorporation.tsv'), sep = '\t', row.name = FALSE)

  invisible(NULL)
}
