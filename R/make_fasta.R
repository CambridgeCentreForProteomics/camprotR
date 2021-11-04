#' Make a FASTA using UniProt accessions
#'
#' @description Given a vector of UniProt accessions, this function will:
#' 1. Download the sequences
#' 2. Print the current UniProt release (you should put this in the FASTA file name)
#' 3. (Optional) add cRAP numbers to the FASTA headers for each sequence
#' 4. Save the sequences in a FASTA file
#'
#' @param accessions `character vector`, the UniProt accessions to use
#' @param file `character`, filepath to save the fasta to e.g. `"crap.fasta"`
#' @param is_crap `logical`, Is the output going to be a cRAP database? If `TRUE`
#' cRAP001, cRAP002, etc. is appended to the sequence headers in the FASTA file.
#' Default is `FALSE`
#' @param overwrite `logical`, if the FASTA file already exists should it be
#' overwritten? Default is `FALSE`
#' @param verbose `logical`, should the function send messages to the console?
#' Default is `TRUE`
#'
#' @return Returns a FASTA file saved to disk at the specified file path.
#' @examples
#' # specify some UniProt accessions
#' crap <- get_ccp_crap()
#'
#' \dontrun{
#' make_fasta(crap, "2021-01_cRAP.fasta")
#' }
#'
#' @export
make_fasta <- function(accessions, file, is_crap = FALSE, overwrite = FALSE, verbose = TRUE) {
  # define payload used to query uniprot
  payload <- list(
    query = paste(accessions, collapse = " "),
    from = "ACC+ID",
    to = "ACC",
    format = "fasta"
  )

  # query uniprot and save FASTA to disk
  # security level fix for Ubuntu 20.04
  # see https://msmith.de/2020/10/02/httr-curl-ubuntu-20-04.html
  httr_config <- switch(
    Sys.info()["sysname"],
    "Linux" = httr::config(ssl_cipher_list = "DEFAULT@SECLEVEL=1"),
    httr::config()
  )

  response <- httr::with_config(
    config = httr_config,
    httr::GET(
      url = "https://www.uniprot.org/uploadlists/",
      query = payload,
      config = httr::write_disk(path = file, overwrite = overwrite)
    )
  )

  # basic http error handling
  httr::stop_for_status(response, "query UniProt")

  # print UniProtKB release
  if (verbose) message(paste("Downloading from UniProtKB release:",
                             response$headers$`x-uniprot-release`))

  if (is_crap) {
    crap <- Biostrings::readAAStringSet(filepath = file)
    names(crap) <- sub_crap(names(crap))
    Biostrings::writeXStringSet(crap, filepath = file)
  }
}

#' Insert cRAP numbers into a character vector
#'
#' @description This function takes a character vector where each element is
#' roughly in the form of a UniProt header e.g. `"sp|XXXXXX|YYYY_YYYY Text goes here"`
#' and substitutes a cRAP number in place of the first `|` symbol, e.g.
#' `"sp|cRAP001|XXXXXX|YYYY_YYYY Text goes here"`.
#'
#' @param x `character vector`, each element must have two `|` symbols with
#' some text in between e.g. `|sometext|`.
#' @param start `numeric`, the number to increment from, default is `1`
#' @param width `numeric`, how many digits the cRAP number should be, default is `3`
#'
#' @return Returns a `character vector` the same length as x.
#'
#' @examples
#' # basic use
#' sub_crap(c("|sometext|", "|moretext|"))
#'
#' # start from a different number
#' sub_crap(c("|sometext|", "|moretext|"), start = 88)
#'
#' # increase number width
#' sub_crap(c("|sometext|", "|moretext|"), start = 1111, width = 4)
#'
#' @export
sub_crap <- function(x, start = 1, width = 3) {
  # check input
  if (!all(grepl("\\|.*\\|", x))) {
    stop("This function only works when each element of x is a character of the
         form '|anything|'")
  }

  output <- mapply(
    function(input, num) {sub("\\|", paste0("|cRAP", num, "|"), input)},
    x,
    formatC(seq.int(from = start, to = start + length(x) - 1, by = 1),
            width = width, format = "d", flag = "0"),
    USE.NAMES = FALSE
  )
  output
}

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
  # security level fix for Ubuntu 20.04
  # see https://msmith.de/2020/10/02/httr-curl-ubuntu-20-04.html
  httr_config <- switch(
    Sys.info()["sysname"],
    "Linux" = httr::config(ssl_cipher_list = "DEFAULT@SECLEVEL=1"),
    httr::config()
  )

  response <- httr::with_config(
    config = httr_config,
    httr::GET(
      url = "https://www.uniprot.org/uniprot/P60709"
    )
  )

  # basic http error handling
  httr::stop_for_status(response, "query UniProt")

  # get release
  response$headers$`x-uniprot-release`
}
