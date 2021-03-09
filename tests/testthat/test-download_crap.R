context("download_crap")

#### Setup ---------------------------------------------------------------------

download_crap_helper <- function(type) {
  tmp <- tempfile(fileext = ".fasta")
  on.exit(file.remove(tmp), add = TRUE, after = TRUE)
  download_crap(tmp, type = type)

  # read the FASTA
  crap <- Biostrings::fasta.index(tmp)

  # select columns of interest
  crap[, c("desc", "seqlength")]
}

#### Tests ---------------------------------------------------------------------
test_that("check_uniprot_release works", {
  expect_true(grepl("[0-9]{4}_[0-9]{2}", check_uniprot_release()))
})

test_that('download_crap(type = "ccp") works', {
  expect_equal_to_reference(
    download_crap_helper(type = "ccp"),
    "reference/download_crap.rds"
  )
})
