#' Append sequences to end of cRAP FASTA
#'
#' @description This function is used to add sequences from a FASTA file onto
#' the end of an existing cRAP FASTA file.
#'
#' @param file `character`, file path of FASTA to append
#' @param crap_file `character`, file path of existing cRAP FASTA to append to
#' @param add_crap add_crap `logical`, should cRAP001, cRAP002, etc. be appended to the
#' sequence headers in the FASTA to append? Default is `TRUE`
#' @param add_crap_start `numeric`, what number should the cRAP00x start at?
#' Default is 1.
#'
#' @return Returns the existing cRAP FASTA file with some more sequences added
#' to the end.
#' @examples
#' # Add some commerical protease sequences onto the end of CCP cRAP FASTA
#' \dontrun{
#' append_crap_fasta(
#'   file = "commercial-proteases.fasta",
#'   crap_file = "ccp_crap_2021-01.fasta",
#'   add_crap = TRUE
#' )
#' }
#'
#' @export
append_crap_fasta <- function(file, crap_file, add_crap = TRUE, add_crap_start = 1) {
  # read in file to append
  to_append <- Biostrings::readAAStringSet(file)

  # add cRAP to headers if specified
  if (add_crap) {
    names(to_append) <- sub_crap(names(to_append), start = add_crap_start)
  }

  # append sequences
  Biostrings::writeXStringSet(to_append, filepath = crap_file, append = TRUE)
}
