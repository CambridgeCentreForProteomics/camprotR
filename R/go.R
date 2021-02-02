#' Determine GO ancestor object
#'
#' @description For a single GO term, get its ontology and then return the
#' correct ancestor mapping object for that ontology.
#'
#' @param term `string`. A single GO.ID.
#' @param ontology `string`. A single ontology, one of: BP, MF, CC.
#'
#' @return Returns an `AnnDbBimap` object which maps GO terms (BP, MF, or
#' CC) to all ancestor terms. Each GO term is mapped to a vector of ancestor
#' GO terms.
#' @keywords internal
#' @export
determine_ancestor_function <- function(term, ontology) {
  ontology <- ontology
  if (is.na(ontology)) {
    return(NA)
  }
  if (ontology == "MF") {
    return(GO.db::GOMFANCESTOR)
  }
  else if (ontology == "CC") {
    return(GO.db::GOCCANCESTOR)
  }
  else if (ontology == "BP") {
    return(GO.db::GOBPANCESTOR)
  }
  else {
    return(NA)
  }
}

#' Determine GO offspring object
#'
#' @description For a single GO term, get its ontology and then return the
#' correct offpsring mapping object for that ontology.
#'
#' @param term `string`. A single GO.ID.
#' @param ontology `string`. A single ontology, one of: BP, MF, CC.
#'
#' @return Returns an `AnnDbBimap` object which maps GO terms (BP, MF, or
#' CC) to all offspring terms. Each GO term is mapped to a vector of offsrping
#' GO terms.
#' @keywords internal
#' @export
determine_offspring_function <- function(term, ontology) {
  if (is.na(ontology)) {
    return(NA)
  }
  if (ontology == "MF") {
    return(GO.db::GOMFOFFSPRING)
  }
  else if (ontology == "CC") {
    return(GO.db::GOCCOFFSPRING)
  }
  else if (ontology == "BP") {
    return(GO.db::GOBPOFFSPRING)
  }
  else {
    return(NA)
  }
}

#' Get all mappings for GO terms
#'
#' For a set of GO terms, obtain all of the ancestor or offspring terms.
#'
#' @param go_ids `character vector`. The GO.IDs to use
#' @param ontologies `named character vector`. Names = GO.IDs and values =
#' Ontologies e.g. BP, CC, MF.
#' @param verbose `logical`.
#' @param direction `string` Either `"ancestor"` or
#' `"offspring"`.
#'
#' @return Returns a `named list` of `character vectors`. Names = GO.IDs and
#' values = all ancestor GO.IDs.
#' @importFrom purrr map2
#' @keywords internal
#' @export
get_all_mappings <- function(go_ids, ontologies, verbose = TRUE, direction = "ancestor") {
  go2relations <- list()

  if (direction == "ancestor") {
    determineFunction <- determine_ancestor_function
  }

  else if (direction == "offspring") {
    determineFunction <- determine_offspring_function
  }

  else {
    stop("direction must be `ancestor` or `offspring`")
  }

  if (verbose) {
    print(sprintf("Getting all %s GO terms for %i observed terms. This may take a while!", direction, length(go_ids)))
  }

  result <- map2(
    ontologies,
    go_ids,
    function(y, x) {

      go_relations <- AnnotationDbi::get(x, determineFunction(x, y))

      if (class(go_relations) == "character") {
        go2relations[[x]] <- go_relations
      } else {
        if (verbose) {
          cat(sprintf("Could not get %s GO terms for GO.ID=%s\n", direction, x))
        }
        go2relations[[x]] <- c(NA)
      }
    }
  )
  return(result)
}

#' Expand data.frame GO terms
#'
#' Input is a data.frame with all the GO terms annotated to a single protein in
#' a single column. Output is a new data.frame with all of the GO terms for
#' that protein (annotated and ancestor).
#'
#' @param go_df `data.frame` for a single protein with a column ==
#' `"GO.ID"`.
#' @param go_col `variable`. Name of column from the data.frame that
#' contains the GO terms.
#' @param go2Ancestor `named list`. Returned by get_all_mappings.
#' Names == GO.IDs and values == all ancestor GO.IDs.
#'
#' @return Returns a `data.frame` containing all ancestor GO terms for a
#' single protein.
#' @keywords internal
#' @export
expand_terms <- function(go_df, go_col, go2Ancestor) {
  observed_go_ids <- unique(go_df[[go_col]]) # unique() shouldn't be req. but does not harm

  unprocessed_ids <- observed_go_ids
  all_ancestors <- observed_go_ids

  while (length(unprocessed_ids) > 0) {
    query_ancestors <- unique(go2Ancestor[unprocessed_ids[1]])
    all_ancestors <- c(all_ancestors, unlist(query_ancestors))

    # removed processed go id
    # and any which were found in ancestors (we don't want to re-query these)
    unprocessed_ids <- unprocessed_ids[-1]
    unprocessed_ids <- base::setdiff(unprocessed_ids, query_ancestors)
  }
  all_ancestors <- setdiff(unique(all_ancestors), "all")

  return(data.frame("GOID" = all_ancestors))
}

#' Get all ancestor GO terms
#'
#' Given a data.frame with with a column containing GO terms, this function
#' will output a data.frame with with the ontologies of those GO terms, their
#' description, and all their ancestor terms.
#'
#' @param go_df `data.frame`. Contains all initial GO terms for proteins
#' of interest with `columns == [feature_col, go_col]`.
#' @param feature_col `string`. Name of the column with the features, e.g.
#' `"UNIPROTKB"`.
#' @param go_col `string`. Name of the column with the GO ids, e.g. `"GO.ID"`.
#' @param return_early_debug `logical`. Stop function early and return a
#' list of GO terms and ancestor GO terms for debugging purposes.
#' @param verbose `logical`.
#'
#' @return Returns a `data.frame` containing all ancestor GO terms for all
#' proteins plus GO term descriptions and ontologies.
#' @importFrom stats setNames
#' @importFrom utils head
#' @export
get_ancestor_go <- function(go_df, feature_col = "UNIPROTKB", go_col = "GO.ID",
                          return_early_debug = FALSE, verbose = TRUE) {

  all_observed_go <- unique(go_df[[go_col]])
  all_observed_go <- all_observed_go[!is.na(all_observed_go)]

  ontologies <- AnnotationDbi::select(
    GO.db::GO.db,
    all_observed_go,
    columns = c("ONTOLOGY"),
    keytype = "GOID"
  )

  ontologies <- setNames(ontologies$ONTOLOGY, ontologies$GOID)
  if (any(is.na(ontologies))) {
    # get names of GO terms that failed to be matched
    bad.terms <- names(ontologies[is.na(ontologies) == TRUE])

    warning(
      "Warning: the following GO terms failed to match to a primary ontology ",
      "in GO.db and were removed: ",
      paste(bad.terms, collapse = ", ")
    )
    # remove NA ontologies
    ontologies <- ontologies[!is.na(ontologies)]

    # remove corresponding GO ids
    all_observed_go <- all_observed_go[!all_observed_go %in% bad.terms]
  }

  go2Ancestor <- get_all_mappings(
    all_observed_go,
    ontologies,
    direction = "ancestor",
    verbose = verbose
  )

  if (return_early_debug) {
    return(list("go_df" = go_df, "go2Ancestor" = go2Ancestor))
  }

  if (verbose) message("Expanding GO terms to include all ancestors for all entries")

  full_go_df <- go_df %>%
    filter(!is.na(all_of(go_col))) %>%
    group_by(.data[[feature_col]]) %>%
    summarise(expand_terms(across(), go_col, go2Ancestor), .groups = "drop_last")

  full_go_details <- AnnotationDbi::select(
    GO.db::GO.db,
    as.character(unique(full_go_df$GOID)),
    columns = c("TERM", "ONTOLOGY"),
    keytype = "GOID"
  )

  full_go_df <- merge(full_go_df, full_go_details, by = "GOID", all.x = TRUE)

  full_go_df <- full_go_df[, c(feature_col, "GOID", "TERM", "ONTOLOGY")]
  colnames(full_go_df)[colnames(full_go_df) == "GOID"] <- go_col

  return(full_go_df)
}
