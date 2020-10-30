context("go")

#### Setup ---------------------------------------------------------------------

# RNA binding and tRNA binding GO terms (both biological process)
rbp.go <- c("GO:0003723", "GO:0000049")

# Make a small GO table of the expected format
go.df <- data.frame(
  UNIPROTKB = c("Q9BYJ9", "Q9NWX6"),
  GO.ID = rbp.go
)

#### Tests ---------------------------------------------------------------------

test_that("determine_ancestor_function works", {
  go.term <- rbp.go[1]

  go.relations <- AnnotationDbi::get(go.term, determine_ancestor_function(go.term, "MF"))
  expect_equal(length(go.relations), 6)
  expect_equal(go.relations[6], "all")
})

test_that("determine_offspring_function works", {
  go.term <- rbp.go[1]

  go.relations <- AnnotationDbi::get(go.term, determine_offspring_function(go.term, "MF"))
  expect_equal(length(go.relations), 173)
})

test_that("get_all_mappings works", {
  rbp.ontologies <- c("MF", "MF")
  names(rbp.ontologies) <- rbp.go

  rbp.mappings <- get_all_mappings(
    go_ids = rbp.go,
    ontologies = rbp.ontologies,
    verbose = TRUE,
    direction = "offspring"
  )

  expect_equal(length(rbp.mappings), 2)
  expect_equal(length(rbp.mappings[[1]]), 173)
  expect_equal(length(rbp.mappings[[2]]), 4)
})

test_that("expand_terms works", {
  rbp.expanded <- expand_terms(go_df = go.df, go_col = "GO.ID", go2Ancestor = rbp.mappings)

  expect_equal(rbp.expanded$GOID, unique(c("GO:0003723", rbp.mappings[[1]], rbp.mappings[[2]])))
})

test_that("get_ancestor_go works", {
  ancestor.terms <- get_ancestor_go(go.df, verbose = FALSE)

  expect_equal(nrow(ancestor.terms), 13)
  expect_false(any(is.na(ancestor.terms)))
})
