#' Extract peptide sequence:protein ID mappings
#'
#' @description Given a list of `infiles`, read the files into a list of
#' data.frames, retain the peptide sequence and protein ID columns, then
#' combine the data.frames into one by binding the rows.
#'
#' @param infiles `list`, paths for PD peptideGroups.txt files to be read
#' by \code{\link[utils]{read.delim}}.
#' @param seq_col `string`, column name for peptide sequence e.g. `"Sequence"`.
#' @param prot_col `string`, column name for protein accessions e.g.
#' `"Protein.Accessions"`.
#' @param master_prot_col `string`, column name for master protein accessions
#' e.g. `"Master.Protein.Accessions"`.
#'
#' @return Returns a `data.frame` object which maps peptide sequences to
#' protein and master protein ID(s) for each infile. Each row is a peptide.
#' @keywords internal
read_infiles <- function(infiles,
                         seq_col,
                         prot_col,
                         master_prot_col) {
  seq2prot <- lapply(infiles, function(x) {
    utils::read.delim(x, stringsAsFactors = FALSE) %>%
      dplyr::select(all_of(c(seq_col, prot_col, master_prot_col))) %>%
      mutate(infile = x)
  }
  ) %>%
    do.call(what = 'rbind')

  return(seq2prot)
}

#' Map peptide sequence to protein across multiple files using
#' approximate parsimony
#'
#' @description Generate a single map of peptide sequences to proteins for
#' multiple files (e.g. multiple separate peptideGroups.txt PD files) using
#' an approximate parsimony approach, specifically the least number
#' of proteins required to account for all observed peptide sequences.
#'
#' @param seq2prot `data.frame`, output from \code{\link{read_infiles}}.
#' @param seq_col `string`, column name for peptide sequence e.g. `"Sequence"`.
#' @param prot_col `string`, column name for protein accessions e.g.
#' `"Protein.Accessions"`.
#' @param master_prot_col `string`, column name for master protein accessions
#' e.g. `"Master.Protein.Accessions"`.
#'
#' @return Returns a `data.frame` object which uniquely maps peptide sequence
#' to master protein(s).
#' @keywords internal
get_sequence2protein <- function(seq2prot,
                                 seq_col,
                                 prot_col) {

  # get unique seq->prot mappings
  seq_protein <- seq2prot %>%
    separate_rows(.data[[prot_col]], sep='; ') %>%
    dplyr::select(all_of(c(seq_col, prot_col))) %>%
    distinct()

  protein_to_sequences <- seq_protein %>%
    group_by_at(prot_col) %>%
    summarise(Sequences = list(!!(sym(seq_col))))

  protein_to_sequences_map <- protein_to_sequences$Sequences

  names(protein_to_sequences_map) <- protein_to_sequences[[prot_col]]

  protein_counts <- protein_to_sequences_map %>%
    names() %>%
    sapply(function(x) length(protein_to_sequences_map[[x]]))

  sequences_left_to_account_for <- unique(seq_protein[[seq_col]])
  sequences_accounted_for <- NULL
  n_seq_per_protein_values <- rev(names(table(protein_counts)))

  s2p <- NULL
  for(n in n_seq_per_protein_values){

    if(length(sequences_left_to_account_for) == 0){
      break()
    }

    proteins <- protein_counts[protein_counts == n]

    mini_s2p <- rep(names(proteins), each = n)
    names(mini_s2p) <- unlist(protein_to_sequences_map[names(proteins)], use.names = FALSE)
    mini_s2p <- mini_s2p[!names(mini_s2p) %in% sequences_accounted_for]

    rep_seq <- intersect(names(mini_s2p), sequences_accounted_for)

    if(length(rep_seq > 0)){
      stop(sprintf('Something has gone wrong here, accounting for the same sequence more than once %s',
                   paste(rep_seq, collapse=',')))
    }

    sequences_accounted_for <- union(sequences_accounted_for, names(mini_s2p))
    sequences_left_to_account_for <- setdiff(sequences_left_to_account_for, names(mini_s2p))

    s2p <- c(s2p, mini_s2p)
  }

  new_seq_to_master <-
    `colnames<-`(data.frame(names(s2p), s2p), c(seq_col, 'Protein')) %>%
    group_by(across(seq_col)) %>%
    summarise('Updated.Master.Protein.Accessions' = paste(.data$Protein, collapse='; '))

  new_seq_to_master
}

#' Compare two protein accession strings
#'
#' @description Check if two protein accession strings are the same,
#' allowing for:
#'
#' 1. multiple accessions (default separated by '; ')
#' 1. multiple accessions in any order
#'
#' @param proteins1 `string`, protein acession(s), set 1
#' @param proteins2 `string`, protein accession(s), set 2
#' @param sep1 `string`, delimiter for proteins1. Default is `'; '`.
#' @param sep2 `string`, delimiter for proteins2. Default is `'; '`.
#'
#' @return Returns `TRUE` if the proteins1 and proteins2 contain the same set
#' of protein accessions.
#' @keywords internal
compare_proteins <- function(proteins1, proteins2, sep1='; ', sep2='; ') {
  if(proteins1 == proteins2){
    return(TRUE)
  }

  proteins1 <- strsplit(proteins1, sep1)[[1]]
  proteins2 <- strsplit(proteins2, sep2)[[1]]

  return(setequal(proteins1, proteins2))
}

#' Plot the number of peptides per protein pre/post-consensus master protein
#' assignment
#'
#' @param seq2new_master_prot `data.frame`, output of merging output of
#' \code{\link{read_infiles}} & \code{\link{get_sequence2protein}}.
#' @param seq_col `string`, column name for peptide sequence e.g. `"Sequence"`.
#' @param master_prot_col `string`, column name for master protein accessions
#' e.g. `"Master.Protein.Accessions"`.
#'
#' @return Returns a `ggplot` showing the number of sequences per master protein
#' in the original infiles, and after the consensus master protein assignment.
#' @keywords internal
plot_sequences_per_protein <- function(seq2new_master_prot,
                                       seq_col,
                                       master_prot_col) {

  master_protein_counts <- seq2new_master_prot %>%
    filter(!!sym(master_prot_col) != '') %>%
    dplyr::select(one_of(c(seq_col, master_prot_col))) %>%
    distinct() %>%
    group_by(across(master_prot_col)) %>%
    tally() %>%
    group_by(n) %>%
    tally() %>%
    mutate(type = 'Original')

  new_master_protein_counts <- seq2new_master_prot %>%
    filter(.data$Updated.Master.Protein.Accessions != '') %>%
    dplyr::select(one_of(c(seq_col, 'Updated.Master.Protein.Accessions'))) %>%
    distinct() %>%
    group_by(.data$Updated.Master.Protein.Accessions) %>%
    tally() %>%
    group_by(n) %>%
    tally() %>%
    mutate(type = 'Updated')

  p <- new_master_protein_counts %>%
    rbind(master_protein_counts) %>%
    arrange(.data$type) %>%
    mutate(type = factor(.data$type, levels = c('Original', 'Updated'))) %>%
    group_by(.data$type) %>%
    mutate(cum_nn = cumsum(.data$nn) / sum(.data$nn)) %>%
    ggplot(aes(x = log10(n), y = .data$cum_nn, colour = .data$type)) +
    geom_step() +
    theme_bw() +
    theme(aspect.ratio = 1) +
    scale_colour_discrete(name = 'Master protein assignment') +
    ylim(0, NA) +
    labs(
      x = 'Sequences per protein (log10)',
      y = 'Fraction of proteins'
    )

  invisible(p)
}

#' Identify the proportion of exact matches in protein accessions
#'
#' @description Identify the proportion of exactly matching protein accessions
#' in a data.frame containing original and updated master protein assignments.
#'
#' @param seq2new_master_prot `data.frame`, output of merging
#' \code{\link{read_infiles}} & \code{\link{get_sequence2protein}}.
#' @param master_prot_col `string`, column name for master protein accessions
#' e.g. `"Master.Protein.Accessions"`.
#'
#' @return Prints the proportion of exact matches in the console.
#' @keywords internal
compare_sequences_per_protein <- function(seq2new_master_prot, master_prot_col){
  id_match <- seq2new_master_prot %>%
    rowwise() %>%
    mutate(same_id=ifelse(compare_proteins(
      !!sym(master_prot_col), .data$Updated.Master.Protein.Accessions),
      'Same ID(s)', 'Different ID(s)')) %>%
    pull(.data$same_id) %>%
    table()

  print(id_match)
  print(round(100*id_match/sum(id_match), 1))

  invisible(NULL)
}



#' Identify the proportion of peptides with a single assigned master protein
#'
#' @description Identify the proportion of peptides with a single assigned
#' master protein in a data.frame containing original and updated master protein
#' assignments.
#'
#' @param seq2new_master_prot `data.frame`, output of merging
#' \code{\link{read_infiles}} & \code{\link{get_sequence2protein}}.
#' @param master_prot_col `string`, column name for master protein accessions
#' e.g. `"Master.Protein.Accessions"`.
#'
#' @return Returns a `ggplot` showing the number of peptide sequences per
#' master protein  in the original infiles, and after the consensus master
#' protein assignment.
#' @keywords internal
compare_single_master <- function(seq2new_master_prot,
                                  master_prot_col) {

  single_master <- seq2new_master_prot %>%
    mutate(original = ifelse(grepl('; ', !!(sym(master_prot_col))),
                             'Original:Multiple',
                             'Original:Single'),
           updated = ifelse(grepl('; ', .data$Updated.Master.Protein.Accessions),
                            'Updated:Multiple',
                            'Updated:Single')) %>%
    group_by(.data$original, .data$updated) %>%
    count() %>%
    spread(key = .data$updated, value = n) %>%
    tibble::column_to_rownames('original')

  print(single_master)
  print(round(100*single_master / sum(single_master), 2))

  invisible(NULL)
}

#' Get consensus peptide to master protein mapping across multiple files
#'
#' @description Given multiple peptideGroups.txt files, identify the
#' parsimonious master proteins to account for the all the peptides identified
#' and (optionally) compare the new and old master protein assignments.
#'
#' The output of this function can then be used when merging multiple
#' peptideGroups.txt files (i.e. the outputs of separate PD studies).
#'
#' @param infiles `list`, paths for PD peptideGroups.txt files to be read
#' by \code{\link[utils]{read.delim}}.
#' @param seq_col `string`, column name for peptide sequence e.g. `"Sequence"`.
#' @param prot_col `string`, column name for protein accessions e.g.
#' `"Protein.Accessions"`.
#' @param master_prot_col `string`, column name for master protein accessions
#' e.g. `"Master.Protein.Accessions"`.
#' @param compare_old_new `logical`, compare the consensus master protein
#' assignments with the original assignments. Default is `TRUE`.
#'
#' @return Returns a `data.frame` which maps peptide sequences to
#' consensus master protein accession(s). Each row is a peptide.
#' @examples
#' \dontrun{
#' # The input in this example is 2x peptideGroups.txt outputs from Proteome
#' # Discoverer. Both files must have the columns (cols) specified. This may require
#' # re-exporting the peptideGroups.txt file from PD as the "Protein.Accessions"
#' # column is often not there by default.
#' get_parsimony_pep2prot(
#'   list("peptideGroups1.txt", "peptideGroups2.txt"),
#'   seq_col = "Sequence",
#'   prot_col = "Protein.Accessions",
#'   master_prot_col = "Master.Protein.Accessions",
#'   compare_old_new = TRUE
#' )
#'
#' }
#' @export
get_parsimony_pep2prot <- function(infiles,
                                   seq_col = 'Sequence',
                                   prot_col = 'Protein.Accessions',
                                   master_prot_col = 'Master.Protein.Accessions',
                                   compare_old_new = TRUE) {
  seq2prot <- read_infiles(infiles, seq_col, prot_col, master_prot_col)

  new_seq_to_master <- get_sequence2protein(seq2prot, seq_col, prot_col)

  seq2new_master_prot <- merge(seq2prot, new_seq_to_master, by=seq_col)

  if(compare_old_new) {

    message(
      sprintf(
        'With original assignments, %s master proteins. With update, %s master proteins\n',
        length(unique(seq2new_master_prot[[master_prot_col]])),
        length(unique(seq2new_master_prot$Updated.Master.Protein.Accessions)))
    )

    print(plot_sequences_per_protein(seq2new_master_prot, seq_col, master_prot_col))

    message('Comparing Master Protein IDs\n')
    compare_sequences_per_protein(seq2new_master_prot, master_prot_col)
    compare_single_master(seq2new_master_prot, master_prot_col)
  }

  invisible(new_seq_to_master)
}
