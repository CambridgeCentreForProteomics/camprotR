#' Check the current UniProt release
#'
#' @description UniProt releases are published approximately every 8 weeks.
#' This function checks what the current UniProt release is.
#'
#' @return Returns a `character`, the current release number in the format
#' YYYY_XX where YYYY is the calendar year and XX a 2-digit number that
#' is incremented for each release of a given year, e.g. 2010_01, 2010_02, etc.
#'
#' @examples
#' # print release number to console
#' check_uniprot_release()
#'
#' # save release number and use in e.g. a file name
#' rls <- check_uniprot_release()
#'
#' paste0("folder/filename_", rls, ".fasta")
#'
#' @export
check_uniprot_release <- function() {
  # query H. sapiens actin
  response <- httr::GET(url = "https://www.uniprot.org/uniprot/P60709")

  # basic http error handling
  httr::stop_for_status(response, "query UniProt")

  # get release
  response$headers$`x-uniprot-release`
}

#' Download cRAP/contaminant FASTA files
#'
#' @description This function is used to download a specific cRAP/contaminants
#' sequence database to a FASTA file. The different types of cRAP databases
#' available through this function currently include:
#'
#' - **ccp** = the contaminants database used by the Cambridge Centre for
#' Proteomics. It is largely based off of The Global Proteome Machine's
#' cRAP database found [here](https://www.thegpm.org/crap/), with some extra
#' sequences including commercial proteases added to the end.
#'
#' @param file `character`, filepath to save the fasta to e.g. `"crap.fasta"`
#' @param type `character`, type of cRAP database to download e.g. `"ccp"`
#' @param add_crap `logical`, should cRAP001, cRAP002, etc. be appended to the
#' sequence headers in the FASTA file? Default is `TRUE`
#' @param overwrite `logical`, if the FASTA file already exists should it be
#' overwritten? Default is `FALSE`
#' @param verbose `logical`, should the function send messages to the console?
#' Default is `TRUE`
#'
#' @return Returns a file saved to disk at the specified filepath
#' @examples
#' \dontrun{
#' download_crap("path/to/file/crap_2021_01.fasta", type = "ccp")
#' }
#'
#' @export
download_crap <- function(file, type = "ccp", add_crap = TRUE,
                          overwrite = FALSE, verbose = TRUE) {
  if (type == "ccp") {
    download_ccp_crap(file = file, overwrite = overwrite, verbose = verbose)
  } else {
    stop('Only type = "ccp" can be used for now.')
  }

  if (add_crap) {
    crap <- Biostrings::readAAStringSet(filepath = file)
    old_names <- names(crap)
    new_names <- mapply(
      function(on, num) {sub("\\|", paste0("|cRAP", num, "|"), on)},
      old_names,
      formatC(seq.int(old_names), width = 3, format = "d", flag = "0"),
      USE.NAMES = FALSE
    )
    names(crap) <- new_names
    Biostrings::writeXStringSet(crap, filepath = file)
  }
}

#' @noRd
download_ccp_crap <- function(file, overwrite = FALSE, verbose = TRUE) {
  # load ccp crap and extract accessions
  input <- Biostrings::readAAStringSet(
    filepath = system.file("extdata", "cRAP_20190401.fasta", package = "camprotR")
  ) %>%
    names()

  accessions <- regmatches(
    input,
    gregexpr("(?<=\\|)[A-Z,0-9]{6}(?=\\|)", input, perl = TRUE)
  ) %>%
    unlist()

  # define payload used to query uniprot
  payload <- list(
    query = paste(accessions, collapse = " "),
    from = "ACC+ID",
    to = "ACC",
    format = "fasta"
  )

  # query uniprot and save FASTA to disk
  response <- httr::GET(
    url = "https://www.uniprot.org/uploadlists/",
    query = payload,
    config = httr::write_disk(path = file, overwrite = overwrite)
  )

  # basic http error handling
  httr::stop_for_status(response, "query UniProt")

  # print UniProtKB release
  if (verbose) message(paste("Downloading from UniProtKB release:",
                             response$headers$`x-uniprot-release`))

  # append commercial sequences
  # Endoproteinase GluC = UPI000C217376
  # recombinant Lys-C (rLys-C) = can't find sequence in UniParc
  com_seqs <- Biostrings::readAAStringSet(
    filepath = system.file("extdata", "commercial_reagents.fasta", package = "camprotR")
  )
  Biostrings::writeXStringSet(com_seqs, filepath = file, append = TRUE)
}
