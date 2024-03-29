% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ptm.R
\name{add_PTM_positions}
\alias{add_PTM_positions}
\title{Add a column describing the position(s) of the PTM(s) with respect to the protein}
\usage{
add_PTM_positions(
  obj,
  proteome_fasta,
  master_protein_col = "Master.Protein.Accessions",
  sequence_col = "Sequence"
)
}
\arguments{
\item{obj}{\code{data.frame} with PD output at PSM/peptide level}

\item{proteome_fasta}{\code{string} Filepath for proteome fasta}

\item{master_protein_col}{\code{string} Name of column containing master protein IDs}

\item{sequence_col}{\code{string}. Name of column containing peptide sequences}
}
\value{
\code{data.frame}
}
\description{
Identify the position(s) of the PTM(s) with respect to the protein.
This is acheieved by finding the position of the peptide sequence in the protein
and using the position(s) of the PTM(s) in the peptide sequence. Where a sequence
has multiple master proteins or the position(s) of the PTM(s) are unknown, the
position(s) of the PTM(s) with respect to the protein is undefined (NA). Input
is PD output at PSM/peptide level having been passed through \code{parse_PTM_scores}
to add 'filtered_pos' column
}
