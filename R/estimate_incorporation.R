#' Estimate SILAC (2-plex) incorporation rate from PD output
#'
#' @description SILAC (2-plex) incorporation can be estimated using the PSM and
#' peptide groups PD outputs, and summarised at the peptide or protein level.
#'
#' This function takes the PSM.txt and PeptideGroups.txt files as inputs and
#' outputs 3 plots and a .tsv file into the designated output directory.
#'
#' @details **Peptide sequencing/mass shift identification**
#' In SILAC, the peptide identity of ions can be established through MS2
#' fragmentation ('peptide sequencing') or by mass shift (within tolerance
#' limits) relative to a sequenced peptide.
#'
#' For 2-plex SILAC experiments, ideally, the correlation between Heavy and
#' Light peptide intensities is the same regardless of whether Light or Heavy
#' peptides are identified by mass shift or sequencing. If e.g Light peptide
#' intensities are not well correlated when they are identified by mass shift,
#' this may indicate that the mass shift identifications are erroneously
#' picking up ghost peptides', which will make incorporation estimation difficult.
#'
#' **Mixing Heavy and Light material for incorporation rate testing**
#' To get around the issue of 'ghost peptides', one can spike in
#' Light material (at cell or protein extract-level) to the Heavy material
#' being analysed.
#'
#' @param psm_input `string` or `data.frame`. File path to PSM-level PD output or
#' the PSM-level PD output data.frame created using \link[utils]{read.delim}.
#' @param peptide_input `string` or `data.frame`. File path to peptide-level PD output or
#' the peptide-level PD output data.frame created using \link[utils]{read.delim}.
#' @param crap_fasta `string`. File path to cRAP fasta used for PD search
#' @param master_protein_col `string`. Name of column containing master
#' protein accessions.
#' @param protein_col `string`. Name of column containing all protein
#' accessions.
#' @param sequence_col `string`. Name of column containing peptide sequences.
#' @param modifications_col `string`. Name of column containing peptide modifications.
#' @param abundance_col_L `string`. Name of column containing light peptide
#' intensity data, or a regex matching the names of multiple abundance columns.
#' @param abundance_col_H `string`. Name of column containing heavy peptide
#' intensity data, or a regex matching the names of multiple abundance columns.
#' @param mix `numeric`. If Light material has been spiked in,
#' what is the abundance relative to the Heavy material? Default is `mix = 0`
#' e.g. no Light spike in. If they are equal, `mix = 1`.
#' @param outdir `NULL` or `string`. If not `NULL`, file path for directory to
#' save the plots and summary table.
#' @return By default returns a `list` with 3 ggplots (HL correlation, peptide-level
#' incorporation, protein-level incorporation) and 1 summary table. If `outdir`
#' is not `NULL` then the plots and table will be saved into `outdir`.
#' @export
#'
#' @examples
#'
#' \dontrun{
#' # input as file paths
#' estimate_incorporation(
#'   psm_input = "data-raw/Molm_13_P4_PSMs.txt",
#'   peptide_input = "data-raw/Molm_13_P4_PeptideGroups.txt",
#'   crap_fasta = "inst/extdata/cRAP_20190401.fasta.gz",
#'   mix = 1,
#'   outdir = "Molm_13_incorporation/"
#' )
#'
#' # input as data.frames
#' estimate_incorporation(
#'   psm_input = read.delim("data-raw/Molm_13_P4_PSMs.txt"),
#'   peptide_input = read.delim("data-raw/Molm_13_P4_PeptideGroups.txt"),
#'   crap_fasta = "inst/extdata/cRAP_20190401.fasta.gz",
#'   mix = 1,
#'   outdir = "Molm_13_incorporation/"
#' )
#' }
#'
estimate_incorporation <- function(
  psm_input,
  peptide_input,
  crap_fasta,
  master_protein_col = "Master.Protein.Accessions",
  protein_col = "Protein.Accessions",
  sequence_col = "Sequence",
  modifications_col = "Modifications",
  abundance_col_L = "^Abundances.Grouped.(F\\d*.)?Light$",
  abundance_col_H = "^Abundances.Grouped.(F\\d*.)?Heavy$",
  mix = 0,
  outdir = NULL
) {
  # create output directory if it does not already exist
  if (!is.null(outdir)) {if (!dir.exists(outdir)) dir.create(outdir)}

  # type checking of inputs
  stopifnot(class(psm_input) %in% c("data.frame", "character"))
  if (class(psm_input) %in% c("character")) psm_input <- utils::read.delim(psm_input)

  stopifnot(class(peptide_input) %in% c("data.frame", "character"))
  if (class(peptide_input) %in% c("character")) peptide_input <- utils::read.delim(peptide_input)

  # throw an error if peptide input does not contain the required columns
  pep_cols <- c(master_protein_col, protein_col, sequence_col, modifications_col,
                "Quan.Info", "Number.of.Missed.Cleavages")

  if (!all(pep_cols %in% colnames(peptide_input))) {
    stop(
      paste("The peptideGroups input is missing the following required columns:",
            pep_cols[!pep_cols %in% colnames(peptide_input)])
    )
  }

  if (!any(grepl(abundance_col_L, colnames(peptide_input)))) {
    stop("The provided abundance_col_L column(s) are missing or misspelled.")
  }
  if (!any(grepl(abundance_col_H, colnames(peptide_input)))) {
    stop("The provided abundance_col_H column(s) are missing or misspelled.")
  }

  # for each peptide, check whether it was MS2 sequenced and what the maximum
  # isolation interference (%) was across all PSMs for that peptide
  psm_sequenced_data <- silac_psm_seq_int(
    psm_input,
    sequence_col = sequence_col,
    mod_col = modifications_col,
    include_interference = TRUE,
    interference_col = "Isolation.Interference.in.Percent"
  )

  # load the cRAP FASTA used for the PD search
  crap_fasta <- Biostrings::fasta.index(crap_fasta, seqtype = "AA")

  # extract the UniProt accessions associated with each cRAP protein
  crap_accessions <- regmatches(
    crap_fasta$desc,
    gregexpr("(?<=\\|).*?(?=\\|)", crap_fasta$desc, perl = TRUE)
  ) %>%
    unlist()

  # parse peptideGroups.txt and filter out:
  # filter out: cRAP, non-tryptic peptides, peptides with redundant Quan
  peptide_data <- parse_features(
    peptide_input,
    master_protein_col = master_protein_col,
    protein_col = protein_col,
    unique_master = FALSE,
    silac = TRUE,
    level = "peptide",
    filter_crap = TRUE,
    crap_proteins = crap_accessions,
    filter_associated_crap = TRUE
  ) %>%
    filter(.data$Quan.Info != "Redundant",
           grepl("K|R", .data[[sequence_col]], ignore.case = TRUE))

  # rename Heavy and Light abundance columns
  colnames(peptide_data)[
    grepl(abundance_col_L, colnames(peptide_data))
  ] <- "Light"
  colnames(peptide_data)[
    grepl(abundance_col_H, colnames(peptide_data))
  ] <- "Heavy"

  # define the columns we need from peptide data
  peptide_data_cols <- c(
    "Sequence",
    "Modifications",
    "Protein.Accessions",
    "Master.Protein.Accessions",
    "Number.of.Missed.Cleavages",
    "Light",
    "Heavy"
  )

  # subset peptide data
  peptide_data <- subset(peptide_data, select = peptide_data_cols)

  # correct the Light channel based on the mixing ratio
  peptide_data$Light.corrected <- peptide_data$Light - (mix * peptide_data$Heavy)

  # calculate incorporation
  peptide_data$Incorporation <- mapply(
    get_incorporation, peptide_data$Light, peptide_data$Heavy
  )

  peptide_data$Incorporation.corrected <- mapply(
    get_incorporation, peptide_data$Light.corrected, peptide_data$Heavy
  )

  # remove SILAC modifications from the Modifications column of peptide data
  peptide_data$Modifications <- remove_silac_modifications(peptide_data$Modifications, level = "peptide")

  # merge the peptide and PSM data
  merged_data <- merge(peptide_data, psm_sequenced_data, by = c(sequence_col, modifications_col))

  # make columns indicating if peptide was identified by PSM or by mass shift
  merged_data$Found.Light <- ifelse(
    merged_data$matched_Light, "Light: Spectrum matched",
    ifelse(
      is.finite(merged_data$Light), "Light: Detected by mass shift",
      "Light: Not found"
    )
  )

  merged_data$Found.Heavy <- ifelse(
    merged_data$matched_Heavy, "Heavy: Spectrum matched",
    ifelse(
      is.finite(merged_data$Heavy), "Heavy: Detected by mass shift",
      "Heavy: Not found"
    )
  )

  # plot Light vs Heavy intensities for each peptide
  p1 <- merged_data %>%
    filter(is.finite(.data$Light) & is.finite(.data$Heavy)) %>%
    ggplot(aes(x = log10(.data$Heavy), y = log10(.data$Light))) +
    geom_point(size = 0.2) +
    theme_bw() +
    theme(aspect.ratio = 1) +
    facet_wrap(.data$Found.Heavy ~ .data$Found.Light) +
    labs(x = "Heavy intensity (log10)", y = "Light intensity (log10)") +
    coord_cartesian(xlim = c(3, NA), ylim = c(3, NA))

  if (mix > 0) {
    p1 <- p1 +
      geom_abline(
        linetype = 2,
        slope = 1,
        intercept = log10((mix / 2) / (1 - (mix / 2))),
        colour = get_cat_palette(1)
      )
  }

  # plot peptide-level incorporation histogram
  if (mix > 0) {
    p2 <- merged_data %>%
      filter(is.finite(.data$Light) & is.finite(.data$Heavy)) %>%
      plot_incorporation(level = "peptide", mix = mix)
  } else {
    p2 <- plot_incorporation(merged_data, level = "peptide", mix = 0)
  }

  # count the number of unique peptides per master protein
  unique_peptides <- merged_data %>%
    select(tidyselect::all_of(c(master_protein_col, sequence_col))) %>%
    unique()

  # get proteins with min. 2 unique peptides
  n_unique_peptides <- unique_peptides[, master_protein_col] %>%
    table() %>%
    as.data.frame() %>%
    `colnames<-`(c(master_protein_col, "n")) %>%
    subset(n > 1)

  # todo: use MsnSets here and make other combine methods available?

  # group data by master protein
  if (mix > 0) {
    merged_data <- merged_data %>%
      filter(is.finite(.data$Light) & is.finite(.data$Heavy))
  }

  protein_data <- merged_data %>%
    group_by(.data[[master_protein_col]]) %>%
    summarise(across(.cols = c("Incorporation", "Incorporation.corrected"), .fns = stats::median)) %>%
    merge(n_unique_peptides, by = master_protein_col)

  # plot protein-level incorporation histogram
  p3 <- plot_incorporation(protein_data, level = 'protein', mix = mix)

  # create table of peptide- and protein-level data
  t1 <- data.frame(matrix(
    unlist(c(
      p2$incorporation_estimates,
      p3$incorporation_estimates
    )),
    ncol = length(p2$incorporation_estimates), byrow = TRUE
  )) %>%
    setNames(names(p2$incorporation_estimates))

  t1$level <- c("peptide", "protein")

  if (!is.null(outdir)) {
    ggsave(file.path(outdir, "HL_correlation.png"), plot = p1)
    ggsave(file.path(outdir, "peptide_incorporation.png"), plot = p2$p)
    ggsave(file.path(outdir, 'protein_incorporation.png'), p3$p)
    utils::write.table(t1, file.path(outdir, 'incorporation.tsv'), sep = '\t', row.name = FALSE)
  } else {
    out <- list(
      'HL_correlation' = p1,
      'peptide_incorporation' = p2$p,
      'protein_incorporation' = p3$p,
      'incorporation_table' = t1
    )
    out
  }
}
