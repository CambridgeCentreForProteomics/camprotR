context("go")

#### Setup ---------------------------------------------------------------------
# Define the ontologies
onts <- c("BP", "CC", "MF")

# Define some GO terms
go_terms <- c(
  "GO:0006378", # BP, mRNA polyadenylation
  "GO:0035770", # CC, Ribonucleoprotein granule
  "GO:0030246" # MF, Carbohydrate binding
)

rbp_terms <- c(
  "GO:0003723", # MF, RNA binding
  "GO:0000049"  # MF, tRNA binding
)
rbp_ontologies <- c("MF", "MF")
names(rbp_ontologies) <- rbp_terms

# Make a small GO table of the expected format
rbp_df <- data.frame(
  UNIPROTKB = c("Q9BYJ9", "Q9NWX6"),
  GO.ID = rbp_terms
)

#### Tests ---------------------------------------------------------------------

test_that("determine_ancestor_function works", {
  go_relations <- mapply(function(gt, ont) {
    AnnotationDbi::get(
      gt, determine_ancestor_function(gt, ont)
    )
  },
  go_terms, onts)

  expect_named(go_relations, go_terms)
  expect_equivalent(lengths(go_relations), c(23, 9, 3))
})

test_that("determine_offspring_function works", {
  go_relations <- mapply(function(gt, ont) {
    AnnotationDbi::get(
      gt, determine_offspring_function(gt, ont)
    )
  },
  go_terms, onts)

  expect_named(go_relations, go_terms)
  expect_equivalent(lengths(go_relations), c(7, 14, 47))
})

test_that("get_all_mappings works with ancestor direction", {
  rbp_mappings <- get_all_mappings(
    go_ids = rbp_terms,
    ontologies = rbp_ontologies,
    verbose = TRUE,
    direction = "ancestor"
  )
  expect_named(rbp_mappings, rbp_terms)
  expect_equivalent(lengths(rbp_mappings), c(6, 7))
})

test_that("get_all_mappings works with offspring direction", {
  rbp_mappings <- get_all_mappings(
    go_ids = rbp_terms,
    ontologies = rbp_ontologies,
    verbose = TRUE,
    direction = "offspring"
  )

  expect_named(rbp_mappings, rbp_terms)
  expect_equivalent(lengths(rbp_mappings), c(173, 4))
})

test_that("expand_terms works", {
  rbp_mappings <- get_all_mappings(
    go_ids = rbp_terms,
    ontologies = rbp_ontologies,
    verbose = TRUE,
    direction = "offspring"
  )

  rbp_expanded <- expand_terms(
    go_df = rbp_df,
    go_col = "GO.ID",
    go2Ancestor = rbp_mappings
  )

  expect_equal(
    rbp_expanded$GOID,
    unique(c("GO:0003723", rbp_mappings[[1]], rbp_mappings[[2]]))
  )
})

test_that("get_ancestor_go works", {
  ancestor_terms <- get_ancestor_go(rbp_df, verbose = FALSE)

  expect_equal(nrow(ancestor_terms), 13)
  expect_false(any(is.na(ancestor_terms)))
})

test_that("get_ancestor_go early debug works", {
  ancestor_terms <- get_ancestor_go(rbp_df, verbose = FALSE, return_early_debug = TRUE)

  expect_identical(ancestor_terms[[1]], rbp_df)
  expect_identical(
    ancestor_terms[[2]],
    get_all_mappings(
      go_ids = rbp_terms,
      ontologies = rbp_ontologies,
      verbose = TRUE,
      direction = "ancestor"
    )
  )
})

#### Sanity checks -------------------------------------------------------------

test_that("determine_ancestor_function produces NA if ontology is NA", {
  result <- determine_ancestor_function(rbp_terms[1], NA)
  expect_identical(result, NA)
})

test_that("determine_ancestor_function produces NA if ontology is nonsense", {
  result <- determine_ancestor_function(rbp_terms[1], "banana")
  expect_identical(result, NA)
})

test_that("determine_offspring_function produces NA if ontology is NA", {
  result <- determine_offspring_function(rbp_terms[1], NA)
  expect_identical(result, NA)
})

test_that("determine_offspring_function produces NA if ontology is nonsense", {
  result <- determine_offspring_function(rbp_terms[1], "banana")
  expect_identical(result, NA)
})

test_that("get_all_mappings throws an error if the direction is nonsense", {
  expect_error(
    get_all_mappings(
      go_ids = rbp_terms,
      ontologies = rbp_ontologies,
      verbose = TRUE,
      direction = "banana"
    ),
    "direction must be `ancestor` or `offspring`"
  )
})

test_that("get_all_mappings complains if it can't find offspring term, and outputs NA", {
  expect_output(
    eisosome_mappings <- get_all_mappings(
      go_ids = "GO:0036286", # eisosome filament
      ontologies = "CC",
      verbose = TRUE,
      direction = "offspring"
    ),
    "Could not get offspring GO terms for GO.ID=GO:0036286"
  )

  expect_identical(eisosome_mappings[[1]], NA)
})

test_that("get_ancestor_go gives a warning for proteins with invalid GO terms", {
  go_df <- data.frame(
    UNIPROTKB = c(rbp_df$UNIPROTKB, "ProteinX"),
    GO.ID = c(rbp_df$GO.ID, "GO:XXXXXXX")
  )

  expect_warning(
    ancestor_terms <- get_ancestor_go(go_df, verbose = TRUE),
    paste(
      "Warning: the following GO terms failed to match to a primary ontology",
      "in GO.db and were removed: GO:XXXXXXX"
    )
  )
  expect_equal(nrow(ancestor_terms), 14)
  expect_true(any(is.na(ancestor_terms$TERM)))
})
