#' Match concatenated IDs in data.frame
#'
#' Match IDs in a data.frame to a reference table and then take other columns
#' from that reference table and bind them to the original data.frame. Mainly
#' used for matching metadata to UniProt accessions, but will work for any
#' type of string ID e.g. Ensembl IDs.
#'
#' @param data `data.frame`. Has the column to be matched.
#' @param to_match `variable`. Name of column in
#' `data` with the IDs to be matched.
#' @param ref `data.frame`. A reference data.frame with the IDs to use for
#' matching and the new columns to add to `data`.
#' @param match `string`. Name of column in `ref` to use for
#' matching.
#' @param new `character vector`. Name of column(s) in `ref` to bind
#' to `data`.
#' @param regex `string`. Regular expression to use for extracting the IDs
#' from the `to_match` column in `data`.
#' @param collapse `string`. If there are multiple IDs in the
#' \code{to_match} column, how should they be collapsed.
#'
#' @return Returns a transformed `data.frame` or `tibble` (depends on structure
#' of input data).
#' @examples
#' ref_df <- data.frame(
#'   accession = c("AAA111", "BBB222", "CCC333", "DDD444"),
#'   name = c("protein a", "protein b", "protein c", "protein d"),
#'   value = c(11, 22, 33, 44)
#' )
#'
#' my_df <- data.frame(
#'   uniprot.id = c("AAA111", "CCC333;BBB222", "EEE555"),
#'   r1 = c(1, 23, 5),
#'   r2 = c(1, 23, 5),
#'   r3 = c(1, 23, 5)
#' )
#'
#' my_df2 <- match_id(
#'   my_df,
#'   uniprot.id,
#'   ref_df,
#'   "accession",
#'   c("name", "value")
#' )
#'
#' @export
match_id <- function(data, to_match, ref, match, new,
                     regex = "[^;]+", collapse = ";") {
  match_col <- eval(substitute(to_match), envir = data)
  result <- match_id_(
    match_col, ref = ref, match = match, new = new,
    regex = regex, collapse = collapse, simplify = FALSE
  )
  do.call(cbind, list(data, result))
}

#' Match concatenated IDs in vector
#'
#' Match IDs in a vector to a reference table and then take other columns
#' from that reference table and output them as a list of vectors. Mainly
#' used for matching metadata to UniProt accessions, but will work for any
#' type of string ID e.g. Ensembl IDs.
#'
#' @param to_match `character vector`. Contains the IDs to be matched.
#' @param ref `data.frame`. A reference data.frame with the IDs to use for
#' matching and the new columns to output as a list.
#' @param match `string`. Name of column in `ref` to use for
#' matching.
#' @param new `character vector`. Name of column(s) in `ref` to output.
#' @param regex `string`. Regular expression to use for extracting the IDs
#' from the `to_match` vector.
#' @param collapse `string`. If there are multiple IDs in the
#' \code{to_match} vector, how should they be collapsed.
#' @param simplify `logical`. Should the output list be unlisted? Default is `FALSE`.
#'
#' @return Returns a list of named vectors unless `simplify = TRUE` wherein a
#' named vector is returned.
#'
#' @examples
#' ref_df <- data.frame(
#'   accession = c("AAA111", "BBB222", "CCC333", "DDD444"),
#'   name = c("protein a", "protein b", "protein c", "protein d"),
#'   value = c(11, 22, 33, 44)
#' )
#'
#' my_vec <- c("AAA111", "CCC333;BBB222", "EEE555")
#'
#' my_df2 <- match_id_(
#'   my_vec,
#'   ref_df,
#'   "accession",
#'   c("name", "value")
#' )
#'
#' @export
match_id_ <- function(to_match, ref, match, new,
                      regex = "[^;]+", collapse = ";", simplify = FALSE) {
  # extract identifiers based on regex
  match_list <- regmatches(
    to_match,
    gregexpr(regex, to_match)
  )

  # match identifiers to the reference table match column
  matched_list <- lapply(match_list, function(x) match(x, ref[[match]]))

  # extract the new columns we want
  new_cols <- sapply(
    new, function(x) {
      Map("[", list(as.character(ref[[x]])), matched_list)
    }, simplify = FALSE
  )

  # collapse multiple matches so the length is the same as match_col
  result <- lapply(new_cols, function(x) sapply(x, function(y) paste(y, collapse = ";")))

  # convert "NA" to NA
  result <- lapply(
    result, function(x) {
      x[x == "NA"] <- NA
      x
    }
  )

  if (simplify) {
    unlist(result)
  } else {
    result
  }
}
