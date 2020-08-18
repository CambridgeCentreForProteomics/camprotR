#' Parse Proteome Discoverer output
#'
#' @description This function parses the output .txt files (peptide groups or PSMs) from
#' Proteome Discoverer and then filters out features based on various criteria.
#'
#' The function performs the following steps:
#'
#' 1. Remove features without a master protein
#' 2. (Optional) Remove features without a unique master protein  (i.e.
#'    Number.of.Protein.Groups == 1)
#' 3. (Optional) Remove features matching a cRAP protein
#' 4. (Optional) Remove features matching any protein associated with
#'   a cRAP protein (see below)
#' 5. Remove features without quantification values (only if TMT or SILAC
#'   are `TRUE` and `level = "peptide"`.)
#'
#' @details **Associated cRAP proteins** are proteins which have at least
#' one feature shared with a cRAP protein. It has been observed that the cRAP
#' database does not contain all possible cRAP proteins e.g. some features
#' can be assigned to a keratin which is not in the provided cRAP database.
#'
#' Using `filter_associated_crap = TRUE` will filter out f2 and f3 in
#' addition to f1, in the example below; regardless of the value in the
#' Master.Protein.Accession column.
#'
#' ```
#' feature  Protein.Accessions         Master.Protein.Accessions
#' f1       protein1, protein2, cRAP,  protein1,
#' f2       protein1, protein3         protein3,
#' f3       protein2                   protein2
#' ```
#'
#' @param data `data.frame` generated from txt file output from Proteome
#' Discoverer.
#' @param master_protein_col `string`. Name of column containing master
#' proteins.
#' @param protein_col `string`. Name of column containing all protein
#' matches.
#' @param unique_master `logical`. Filter out features without a unique
#' master protein.
#' @param silac `logical`. Is the experiment a SILAC experiment?
#' @param TMT `logical`. Is the experiment a TMT experiment?
#' @param level `string`. Type of input file, must be one of either
#' `"peptide"` or `"PSM"`.
#' @param filter_crap `logical`. Filter out features which match a cRAP
#' protein.
#' @param crap_proteins `character vector`. Contains the cRAP accessions,
#' for example: `c("P02768")` which is serum albumin.
#' @param filter_associated_crap `logical`. Filter out features which
#' match a cRAP associated protein.
#' @return `data.frame` with the filtered Proteome Discoverer output.
#' @example inst/examples/ex-parse_features.R
#' @export
parse_features <- function(data,
                           master_protein_col = "Master.Protein.Accessions",
                           protein_col = "Protein.Accessions",
                           unique_master = TRUE,
                           silac = FALSE,
                           TMT = FALSE,
                           level = "peptide",
                           filter_crap = TRUE,
                           crap_proteins = NULL,
                           filter_associated_crap = TRUE) {
  # check arguments
  stopifnot(level %in% c("PSM", "peptide"))
  if (filter_crap) {
    if (is.null(crap_proteins)) {
      stop("must supply the crap_proteins argument to filter cRAP proteins")
    }
  }

  # print input summary
  message("Parsing features...")
  message_parse(data, master_protein_col, "Input")

  # remove crap proteins
  if (filter_crap) {
    message(sprintf("%s cRAP proteins supplied",
                    length(crap_proteins)))

    # identify associated crap proteins first if necessary
    if (filter_associated_crap) {
      associated_crap <- data %>%
        filter(master_protein_col %in% crap_proteins |
                 grepl("cRAP", data[[protein_col]], ignore.case = FALSE)) %>%
        pull(protein_col) %>%
        strsplit("; ") %>%
        unlist()
      associated_crap <- associated_crap[!grepl("cRAP", associated_crap)]

      message(sprintf("%s proteins identified as 'cRAP associated'",
                      length(associated_crap)))
    }

    # then remove normal crap proteins
    data <- data %>%
      filter(!master_protein_col %in% crap_proteins &
               !grepl("cRAP", data[[protein_col]], ignore.case = FALSE))
    message_parse(data, master_protein_col, "cRAP features removed")

    # then remove associated crap proteins if necessary
    if (filter_associated_crap) {
      if (length(associated_crap) > 0) {
        # remove isoforms
        associated_crap_no_isoform <- unique(sapply(strsplit(associated_crap, "-"), "[[", 1))
        associated_crap_regex <- paste(associated_crap_no_isoform, collapse = "|")
        data <- data[!grepl(associated_crap_regex, data[[protein_col]]), ]

        message_parse(data, master_protein_col, "associated cRAP features removed")
      }
    }
  }

  # remove features without a master protein
  if (any(is.na(data[[master_protein_col]]))) {
    data <- data[!is.na(data[[master_protein_col]]), ]
    message_parse(data, master_protein_col, "features without a master protein removed")
  }

  # remove features with non-unique master proteins
  if (unique_master) {
    data <- data[data[["Number.of.Protein.Groups"]] == 1, ]
    message_parse(data, master_protein_col, "features with non-unique master proteins removed")
  }

  # remove features with quantification warnings if necessary
  if (silac | TMT & level == "peptide") {
    data <- data[is.na(data[["Quan.Info"]]), ]
    message_parse(data, master_protein_col, "features without quantification removed")
  }
  data
}

#' Filter a PSM-level MSnSet to remove low quality PSMs
#'
#' @description Filter PSMs from TMT quantification to remove the following
#'
#' 1. Missing values (NA) for all tags
#' 2. Interference/co-isolation above a set value (default=100, e.g no filtering)
#' 3. Signal:noise ratio below a set value (default=0, e.g no filtering)
#' 
#' @param obj `MSnSet` PSMs
#' @param inter_thresh `numeric` Maximum allowed interference/co-isolation
#' @param sn_thresh `numeric` Minimum allowed signal:noise threshold
#' @param master_protein_col `string`. Name of column containing master
#' proteins.
#' @param inter_col `string` Name of column containing the interference value
#' @param sn_col `string` Name of column containing the signal:noise value
#' 
#' @return `MSnSet` with the filtered PSMs.
#' @export
filter_TMT_PSMs <- function(obj,
                        inter_thresh=100,
                        sn_thresh=0,
                        master_protein_col='Master.Protein.Accessions',
                        inter_col='Isolation.Interference.in.Percent',
                        sn_col='Average.Reporter.SN',
                        verbose=TRUE){
  
  if(verbose) message("Filtering PSMs...")
  
  obj <- obj[rowSums(is.finite(exprs(obj)))>0,]
  if(verbose) message_parse(fData(obj), master_protein_col, "No quant filtering")
  
  obj <- obj[fData(obj)[[inter_col]]<=inter_thresh,]
  if(verbose) message_parse(fData(obj), master_protein_col, "Co-isolation filtering")
  
  obj <- obj[fData(obj)[[sn_col]]>=sn_thresh,]
  if(verbose) message_parse(fData(obj), master_protein_col, "S:N ratio filtering")
  
  obj
}
