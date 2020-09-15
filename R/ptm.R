#' Parse the PTM probabilities and add new columns with PTM information
#'
#' @description Extract PTM information from the PTM probabilities column in
#' the output from PD and add new columns with this information. Also optionally
#' prints a summary of how many features (PSMs/peptides) pass a given probability
#' threshold
#'
#' @param obj `data.frame` with PD output at PSM/peptide level
#' @param threshold `numeric` If any score is below a set threshold, disregard all putative PTM sites
#' @param ptm_col `character` Columm name for PTM probabilities
#' @param prob_split `character` regex to split PTM probabilities
#' @param collapse_delimiter `character` delimiter for multiple values in output columns
#' @param verbose Set TRUE to print log of events to console
#'
#' @return `data.frame`
#' @export
parse_PTM_scores <- function(obj,
                           threshold = 95,
                           ptm_col = "PhosphoRS.Best.Site.Probabilities",
                           prob_split = '; |: ',
                           collapse_delimiter = ";",
                           verbose = TRUE) {
  if (class(obj) != "data.frame") {
    stop("'obj' must be a data.frame")
  }

  message(sprintf(
    "Removed %s Features where the ptm_col value == `Inconclusive data`",
    sum(obj[[ptm_col]] == "Inconclusive data")
  ))
  obj <- obj[obj[[ptm_col]] != "Inconclusive data", ]


  # initiate vectors with empty string
  # where the PTM(s) pass the threshold, these will be updated
  filtered_ptm_desc <- filtered_ptm_res <- filtered_ptm_pos <- filtered_ptm <- filtered_ptm_score <- rep("", nrow(obj))

  split_probabities <- strsplit(obj[[ptm_col]], split = prob_split)

  log <- list(
    "Total Features" = 0,
    "Total detected PTMFeatures" = 0,
    "Features passing filter" = 0,
    "Features failing filter" = 0,
    "BiPTM/multiPTM Features where some sites fail filter" = 0,
    "Total detected sites" = 0,
    "Sites passing filter" = 0,
    "Sites failing filter" = 0,
    "monoPTM passing filter" = 0,
    "biPTM passing filter" = 0,
    "multiPTM passing filter" = 0,
    "Too many isoforms" = 0
  )

  for (i in seq_along(split_probabities)) {
    log[["Total Features"]] <- log[["Total Features"]] + 1

    peptide_ptm_scores <- split_probabities[[i]]
    if (length(peptide_ptm_scores) == 0) {
      next()
    } # no PTM detected
    if (is.na(peptide_ptm_scores[[1]])) {
      next()
    } # no PTM detected

    log[["Total detected PTMFeatures"]] <- log[["Total detected PTMFeatures"]] + 1

    if (peptide_ptm_scores[[1]] == "Too many isoforms") {
      log[["Too many isoforms"]] <- log[["Too many isoforms"]] + 1
      log[["Features failing filter"]] <- log[["Features failing filter"]] + 1
      next()
    }

    log[["Total detected sites"]] <- log[["Total detected sites"]] + length(peptide_ptm_scores) / 2

    scores <- peptide_ptm_scores[seq(2, length(peptide_ptm_scores), 2)]

    # if any score is below threshold, disregard all putative ptm sites
    if (any(as.numeric(scores) < threshold)) {
      log[["Sites failing filter"]] <- log[["Sites failing filter"]] + length(peptide_ptm_scores) / 2
      log[["Features failing filter"]] <- log[["Features failing filter"]] + 1
      if (any(as.numeric(scores) >= threshold)) {
        log[["BiPTM/multiPTM Features where some sites fail filter"]] <- log[[
        "BiPTM/multiPTM Features where some sites fail filter"]] + 1
      }
      # if we want to handle this differently, can implement an alternative approach here
      # and move the rest of the code below into an else clause
      next()
    }

    log[["Sites passing filter"]] <- log[["Sites passing filter"]] + length(peptide_ptm_scores) / 2
    log[["Features passing filter"]] <- log[["Features passing filter"]] + 1

    if (length(scores) == 1) {
      log[["monoPTM passing filter"]] <- log[["monoPTM passing filter"]] + 1
    }
    else if (length(scores) == 2) {
      log[["biPTM passing filter"]] <- log[["biPTM passing filter"]] + 1
    }
    else {
      log[["multiPTM passing filter"]] <- log[["multiPTM passing filter"]] + 1
    }

    ptms <- peptide_ptm_scores[seq(1, length(peptide_ptm_scores), 2)] # extract the PTMs info
    split_ptms <- unlist(strsplit(ptms, split = '\\(|\\)')) # split to remove parantheses
    modifications <- split_ptms[seq(2, length(split_ptms), 2)] # extract modifications, e.g "phospho"
    positions <- split_ptms[seq(1, length(split_ptms), 2)] # extract the positions, e.g "S6"
    residues <- substr(positions, 1, 1) # extract first element, e.g S
    positions <- sub('.', '', positions) # remove first element and leave position, e.g 6

    # paste together the value, separated by option(collapse_delimiter) and update vectors which will become columns
    filtered_ptm_desc[[i]] <- paste(peptide_ptm_scores, collapse = collapse_delimiter)
    filtered_ptm_res[[i]] <- paste(residues, collapse = collapse_delimiter)
    filtered_ptm_pos[[i]] <- paste(positions, collapse = collapse_delimiter)
    filtered_ptm[[i]] <- paste(modifications, collapse = collapse_delimiter)
    filtered_ptm_score[[i]] <- paste(scores, collapse = collapse_delimiter)
  }

  # add columns
  obj['filtered_PTM_desc'] = filtered_ptm_desc
  obj['filtered_res'] = filtered_ptm_res
  obj['filtered_pos'] = filtered_ptm_pos
  obj['filtered_ptm'] = filtered_ptm
  obj['filtered_score'] = filtered_ptm_score

  if (verbose) {
    for (event in names(log)) {
      message(sprintf("%s: %i", event, log[[event]]))
    }
  }

  return(obj)
}


#' Add a column describing the position(s) of the PTM(s) with respect to the protein
#'
#' @description Identify the position(s) of the PTM(s) with respect to the protein.
#' This is acheieved by finding the position of the peptide sequence in the protein
#' and using the position(s) of the PTM(s) in the peptide sequence. Where a sequence
#' has multiple master proteins or the position(s) of the PTM(s) are unknown, the
#' position(s) of the PTM(s) with respect to the protein is undefined (NA). Input
#' is PD output at PSM/peptide level having been passed through `parse_PTM_scores`
#' to add 'filtered_pos' column
#'
#' @param obj `data.frame` with PD output at PSM/peptide level
#' @param proteome_fasta `character` Filepath for proteome fasta
#' @param master_protein_col `character` Column name for master protein
#'
#' @return `data.frame`
#' @export
add_PTM_positions <- function(obj, proteome_fasta, master_protein_col = "Master.Protein.Accessions") {
  proteome <- Biostrings::readAAStringSet(proteome_fasta)
  names(proteome) <- sapply(base::strsplit(names(proteome), split = "\\|"), "[[", 2)

  combine_peptide_ptm_positions <- function(proteome, protein, sequence, filtered_pos) {
    # Given a master protein(s), sequence and positions of PTM AA with respect to peptide sequence (filtered_pos),
    # return the PTM positions with respect to protein sequence

    if (filtered_pos == "") {
      return(NA)
    }
    return_string <- NULL

    if (length(strsplit(protein, '; ')[[1]]) > 1) {
      return("")
    }

    peptide_starts <- stats::start(Biostrings::matchPattern(sequence, proteome[[protein]]))

    for (p_start in peptide_starts) {
      position_string <- NULL
      for (ptm_p in strsplit(filtered_pos, split = ";")[[1]]) {
        position_string[[ptm_p]] <- p_start + as.numeric(ptm_p) - 1
      }
      return_string[[as.character(p_start)]] <- paste0(position_string, collapse = "|")
    }

    return(paste0(return_string, collapse = ";"))
  }

  obj$ptm_position <- apply(
    obj,
    MARGIN = 1, function(x) combine_peptide_ptm_positions(
      proteome, x[[master_protein_col]], x[["Sequence"]], x[["filtered_pos"]]
    )
  )

  return(obj)
}


#' Get the amino acid sequence around a PTM
#'
#' @description Get the amino acid sequence around a PTM. Will return NA if
#' peptide maps to multiple proteins or has multiple PTMs. If padding extends
#' outside the protein AA sequence, padding will be extended with '_'.
#'
#' @param proteome `XStringSet` as generated by `Biostrings::readAAStringSet`
#' @param protein `character` protein ID
#' @param ptm_position `numeric` position of the PTM in the protein
#' @param pad `numeric` Number of amino acids around PTM
#'
#' @return `character` PSM-centered amino acid sequence
#' @export
get_sequence <- function(proteome, protein, ptm_position, pad = 7) {
  ptm_position <- suppressWarnings(as.numeric(as.character(ptm_position)))

  if (is.na(ptm_position)) {
    return(NA)
  }

  if (grepl("; ", ptm_position)) {
    return(NA)
  }

  if (!protein %in% names(proteome)) {
    return(NA)
  }

  protein_length <- length(proteome[[protein]])

  if (ptm_position > protein_length) {
    warning(sprintf(
      "PTM positions is outside protein sequence! Returning NA. %s: [-%s], PTM: %s",
      protein, protein_length, ptm_position
    ))
    return(NA)
  }

  start_pad <- end_pad <- ""

  start <- ptm_position - pad
  if (start < 0) {
    start_pad <- paste0(rep("_", start * -1), collapse = "")
    start <- 0
  }

  end <- ptm_position + pad
  if (end > protein_length) {
    end_pad <- paste0(rep("_", (end - protein_length)), collapse = "")
    end <- protein_length
  }

  mod_position <- pad + 1

  sequence <- as.character(proteome[[protein]][start:end])
  sequence <- paste0(start_pad, sequence, end_pad)

  sequence <- paste(base::substr(sequence, 1, pad),
    tolower(base::substr(sequence, pad + 1, pad + 1)),
    base::substr(sequence, pad + 2, pad + pad + 1),
    sep = ""
  )

  return(sequence)
}

#' Add a column with amino acid sequence around a PTM
#'
#' @description Add a column with amino acid sequence around a PTM. Value will be
#' NA if peptide maps to multiple proteins or has multiple PTMs. If padding extends
#' outside the protein AA sequence, padding will be extended with '_'. The PTM-
#' centered AA sequence is useful to integrate external databases. Input
#' is PD output at PSM/peptide level having been passed through `parse_PTM_scores`
#' to add 'filtered_pos' column
#'
#' @param obj `data.frame` with PD output at PSM/peptide level
#' @param proteome_fasta `character` Filepath for proteome fasta
#' @param master_protein_col `character` Column name for master protein
#'
#' @return `data.frame`
#' @export
add_site_sequence <- function(obj,
                            proteome_fasta,
                            master_protein_col = "Master.Protein.Accessions") {
  proteome <- Biostrings::readAAStringSet(proteome_fasta)
  names(proteome) <- sapply(base::strsplit(names(proteome), split = "\\|"), "[[", 2)

  obj <- obj %>%
    rowwise() %>%
    mutate(site_seq = get_sequence(
      proteome,
      !!sym(master_protein_col),
      .data$ptm_position
    )) %>%
    data.frame()

  return(obj)
}

#' Add a column describing the position(s) of the peptide sequence with respect
#' to the protein
#'
#' @description Identify the position(s) of the peptide sequence in the protein.
#' Where a sequence has multiple master proteins, the peptide position is
#' undefined (NA).
#'
#' @param obj `data.frame` with PD output at PSM/peptide level
#' @param proteome_fasta `character` Filepath for proteome fasta
#' @param master_protein_col `character` Column name for master protein
#'
#' @return `data.frame`
#' @export
add_peptide_positions <- function(obj,
                                proteome_fasta,
                                master_protein_col = "Master.Protein.Accessions") {
  proteome <- Biostrings::readAAStringSet(proteome_fasta)
  names(proteome) <- sapply(base::strsplit(names(proteome), split = "\\|"), "[[", 2)

  combine_peptide_ptm_positions <- function(proteome, protein, sequence) {

    # Given a master protein(s) and AA sequence,
    # return the AA position with respect to protein sequence
    # If more than one position, possible, return NA

    if (length(strsplit(protein, '; ')[[1]]) > 1) {
      return(c(NA, NA))
    }

    if (!protein %in% names(proteome)) {
      return(c(NA, NA))
    }
    peptide_start <- stats::start(Biostrings::matchPattern(sequence, proteome[[protein]]))
    peptide_end <- stats::end(Biostrings::matchPattern(sequence, proteome[[protein]]))

    if (length(peptide_start) != 1) {
      return(c(NA, NA))
    }

    else {
      return(c(peptide_start, peptide_end))
    }
  }

  obj[, c('peptide_start', 'peptide_end')] <- t(apply(
    obj,
    MARGIN = 1, function(x) combine_peptide_ptm_positions(
      proteome, x[[master_protein_col]], x[["Sequence"]]
    )
  ))

  obj$peptide_start <- as.numeric(obj$peptide_start)
  obj$peptide_end <- as.numeric(obj$peptide_end)

  return(obj)
}
