context("ptm")

#### Setup ---------------------------------------------------------------------
# Objects for testing get_sequence()
proteome <- Biostrings::readAAStringSet(
  system.file("extdata", "reference.fasta", package = "camprotR"))

names(proteome) <- sapply(strsplit(names(proteome), split='\\|'), '[[', 2)

# Objects for testing add_site_sequence()
ass_input <- data.frame('protein'=c('L0R819', 'P98196', 'P14618'),
                        'ptm_position'=c('20', '100', '40'))

ass_expected_output <- ass_input
ass_expected_output$site_seq <- c('VLIPTDNsTPHKEDL', 'PVTSGLPlFFVITVT', NA)

fasta <- system.file("extdata", "reference.fasta", package = "camprotR")

# Objects for testing add_peptide_positions()
app_input <- data.frame('protein'=c('L0R819', 'P98196', 'P14618'),
                        'Sequence'=c('GTRPEDSSVLIPTDNSTPHK',
                                     'QLWPTATER',
                                     'LNFSHGTHEYHAETIK'))

app_expected_output <- app_input
app_expected_output$peptide_start <- c(5, 1099, NA)
app_expected_output$peptide_end <- c(24, 1107, NA)

#### Tests ---------------------------------------------------------------------
test_that("get_sequence works for protein with single PTM", {
  expect_equal(get_sequence(proteome, 'L0R819', '8', pad=7), "MPSRGTRpEDSSVLI")
})

test_that("get_sequence can add padding at the start", {
  expect_equal(get_sequence(proteome, 'L0R819', '4', pad=7), "___MPSRgTRPEDS")
})

test_that("get_sequence can add padding at the end", {
  expect_equal(get_sequence(proteome, 'L0R819', '90', pad=7), "YQEIENLdKTKIKK_")
})

test_that("add_site_sequence works", {
  expect_equal(add_site_sequence(ass_input, fasta, master_protein_col='protein'),
               ass_expected_output)
})

test_that("add_peptide_position works", {
  expect_equal(add_peptide_positions(app_input, fasta, master_protein_col = "protein"),
               app_expected_output)
})

#### Sanity checks -------------------------------------------------------------

test_that("get_sequence produces NA for protein not in proteome", {
  expect_equal(get_sequence(proteome, 'A0R819', '1', pad=7), NA)
})

test_that("get_sequence produces NA for protein with multiple PTMs", {
  expect_equal(get_sequence(proteome, 'L0R819', '40; 120', pad=7), NA)
})

test_that("get_sequence warns if the PTM position is outside protein length", {
  expect_warning(get_sequence(proteome, 'L0R819', '200', pad=7))
})

test_that("combine_peptide_ptm_positions produces NA if protein is not unique", {
  expect_equal(camprotR:::combine_peptide_ptm_positions(proteome, "L0R819; A0R819", "S"), c(NA, NA))
})

test_that("combine_peptide_ptm_positions produces NA if peptide position not unique", {
  expect_equal(camprotR:::combine_peptide_ptm_positions(proteome, "L0R819", "SS"), c(NA, NA))
})
