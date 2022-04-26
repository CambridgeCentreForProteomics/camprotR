#' Add sequences to end of a FASTA file
#'
#' @description This function is used to add sequences from a FASTA file (file1)
#' onto the end of another FASTA file (file2).
#'
#' If file2 is a cRAP database, then you can optionally add cRAP numbers to the
#' headers of file1, starting at the desired number e.g. cRAP127.
#'
#' @param file1 `character`, file path of FASTA to append.
#' @param file2 `character`, file path of FASTA to append to.
#' @param is_crap `logical`, should cRAP numbers e.g. cRAP001, cRAP002, etc.
#' be added to sequence headers of file1? Default is `FALSE`.
#' @param crap_start `numeric`, what number should the cRAPxxx start at?
#' Default is 1 which will produce: cRAP001.
#'
#' @return Overwrites FASTA file2 with some more sequences added
#' to the end.
#'
#' @examples
#' # Add some commercial protease sequences onto the end of CCP cRAP FASTA
#' \dontrun{
#' append_fasta(
#'   file1 = "commercial-proteases.fasta",
#'   file2 = "2021-01_CCP_cRAP.fasta",
#'   is_crap = TRUE,
#'   crap_start = 128
#' )
#' }
#'
#' @export
append_fasta <- function(file1, file2, is_crap = FALSE, crap_start = 1) {
  # read in file to append
  to_append <- Biostrings::readAAStringSet(file1)

  # add cRAP to headers if specified
  if (is_crap) {
    names(to_append) <- sub_crap(names(to_append), start = crap_start)
  }

  # append sequences
  Biostrings::writeXStringSet(to_append, filepath = file2, append = TRUE)
}
