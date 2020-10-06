context("match_id")

#### Setup ---------------------------------------------------------------------

my_df <- data.frame(
  accession = c("A1A1A1;B2B2B2", "C3C3C3", "Z9Z9Z9"),
  sample1 = c(1, 1, 1),
  sample2 = c(2, 2, 2),
  sample3 = c(3, 3, 3)
)

my_df2 <- data.frame(
  accession = c("A1A1A1;B2B2B2", "C3C3C3", "Z9Z9Z9"),
  sample4 = c(4, 4, 4)
)

ref_df <- data.frame(
  entry = c("A1A1A1", "B2B2B2", "C3C3C3", "D4D4D4"),
  gene.name = c("abb1", "bcc2", "cdd3", "dee4"),
  protein.name = c("Protein A", "Protein B", "Protein C", "Protein D")
)

#### Tests ---------------------------------------------------------------------

test_that("match_id defaults work", {
  my_df3 <- match_id(my_df, accession, ref_df, "entry", "gene.name")

  expect_equal(
    my_df3$gene.name, c("abb1;bcc2", "cdd3", NA)
  )
})

test_that("match_id can add multiple columns", {
  my_df3 <- match_id(my_df, accession, ref_df, "entry", c("gene.name", "protein.name"))

  expect_equal(
    my_df3$gene.name, c("abb1;bcc2", "cdd3", NA)
  )

  expect_equal(
    my_df3$protein.name, c("Protein A;Protein B", "Protein C", NA)
  )
})

test_that("match_id alternative regex and type coercion warning works", {
  expect_warning(
    my_df3 <- match_id(my_df, accession, my_df2, "accession", "sample4", regex = ".*")
  )
  expect_equal(
    as.numeric(my_df3$sample4), c(4, 4, 4)
  )
})
