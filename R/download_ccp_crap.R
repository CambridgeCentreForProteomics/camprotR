#' Get the UniProt sequences in CCP cRAP
#'
#' @description This function take no inputs and simply outputs a vector of the
#' UniProt accessions for the sequences in the CCP cRAP FASTA.
#'
#' @return Returns a `character vector` of UniProt accessions.
#' @export
get_ccp_crap <- function() {
  # load ccp crap and extract accessions
  input <- Biostrings::readAAStringSet(
    filepath = system.file("extdata", "cRAP_20190401.fasta.gz", package = "camprotR")
  ) %>%
    names()

  accessions <- regmatches(
    input,
    gregexpr("(?<=\\|)[A-Z,0-9]{6}(?=\\|)", input, perl = TRUE)
  ) %>%
    unlist()

  accessions[!accessions %in% "000000"]

  # deal with alpha amylase proteins specifically
  accessions[4:6] <- c("P0DUB6", "P0DTE7", "P0DTE8")
  accessions
}

#' Download CCP cRAP FASTA file
#'
#' @description This function is used to download the Cambridge Centre for
#' Proteomics cRAP FASTA database using sequences from the lastest UniProt release.
#'
#' @param file `character`, filepath to save the fasta to e.g. `"crap.fasta"`
#' @param is_crap `logical`, should cRAP001, cRAP002, etc. be appended to the
#' sequence headers in the FASTA file? Default is `TRUE`
#' @param overwrite `logical`, if the FASTA file already exists should it be
#' overwritten? Default is `FALSE`
#' @param verbose `logical`, should the function send messages to the console?
#' Default is `TRUE`
#'
#' @return Returns a FASTA saved to disk at the specified file path.
#' @examples
#' \dontrun{
#' download_ccp_crap("path/to/file/2021-01_CCP_cRAP.fasta")
#' }
#'
#' @export
download_ccp_crap <- function(file, is_crap = TRUE, overwrite = FALSE, verbose = TRUE) {
  accessions <- get_ccp_crap()

  make_fasta(accessions = accessions,
                  file = file,
                  is_crap = is_crap,
                  overwrite = overwrite,
                  verbose = verbose)

  # Add commercial sequences on the end
  append_fasta(
    file1 = system.file("extdata", "commercial_reagents.fasta", package = "camprotR"),
    file2 = file,
    is_crap = is_crap,
    crap_start = 126
  )
}
