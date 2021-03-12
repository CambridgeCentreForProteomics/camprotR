context("download_ccp_crap")

#### Setup ---------------------------------------------------------------------

download_ccp_crap_helper <- function() {
  tmp <- tempfile(fileext = ".fasta")
  on.exit(file.remove(tmp), add = TRUE, after = TRUE)
  download_ccp_crap(tmp)

  # read the FASTA
  crap <- Biostrings::fasta.index(tmp)

  # select columns of interest
  crap[, c("desc", "seqlength")]
}

#### Tests ---------------------------------------------------------------------
test_that("check_uniprot_release works", {
  expect_true(grepl("[0-9]{4}_[0-9]{2}", check_uniprot_release()))
})

test_that('download_ccp_crap works', {
  expect_equal_to_reference(
    download_ccp_crap_helper(),
    "reference/download_ccp_crap.rds"
  )
})

test_that("sub_crap works", {
  expect_equal(sub_crap("|text|"), "|cRAP001|text|")
  expect_equal(sub_crap("|text|", 77), "|cRAP077|text|")
  expect_equal(sub_crap("|text|", 111, 4), "|cRAP0111|text|")
})
