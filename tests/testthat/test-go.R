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

test_that("get_ancestor_go works", {
  ancestor.terms <- get_ancestor_go(go.df, verbose = FALSE)

  expect_equal(nrow(ancestor.terms), 13)
  expect_false(any(is.na(ancestor.terms)))
})
