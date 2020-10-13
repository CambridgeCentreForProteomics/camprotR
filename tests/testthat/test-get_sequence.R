context('Get site sequence')

# set up #
proteome <- Biostrings::readAAStringSet(
  system.file("extdata", "reference.fasta", package = "camprotR"))

names(proteome) <- sapply(strsplit(names(proteome), split='\\|'), '[[', 2)
##########

# tests #
test_that("single protein", {
  expect_equal(get_sequence(proteome, 'L0R819', '8', pad=7), "MPSRGTRpEDSSVLI")
})

test_that("single protein, pad outside start", {
  expect_equal(get_sequence(proteome, 'L0R819', '4', pad=7), "___MPSRgTRPEDS")
})

test_that("single protein, pad outside end", {
  expect_equal(get_sequence(proteome, 'L0R819', '90', pad=7), "YQEIENLdKTKIKK_")
})

test_that("multiple positions", {
  expect_equal(get_sequence(proteome, 'A0R819', '1', pad=7), NA)
})

test_that("Protein not in proteome", {
  expect_equal(get_sequence(proteome, 'L0R819', '40; 120', pad=7), NA)
})

test_that("position outside boundary", {
  expect_warning(get_sequence(proteome, 'L0R819', '200', pad=7))
})

