context("ptm")

#### Setup ---------------------------------------------------------------------
# Take 10% of psm_tmt_phospho dataset for testing purposes
phospho <- psm_tmt_phospho[1:500, ]

# Create an objects for testing add_PTM_positions()
nucleolin_psm <- psm_tmt_phospho[1, ]

# Objects for testing get_sequence()
proteome <- Biostrings::readAAStringSet(
  system.file("testdata", "reference.fasta",
              package = "camprotR")
)

names(proteome) <- sapply(strsplit(names(proteome), split='\\|'), '[[', 2)

# Objects for testing add_site_sequence()
ass_input <- data.frame('protein'=c('L0R819', 'P98196', 'P14618'),
                        'ptm_position'=c('20', '100', '40'))

ass_expected_output <- ass_input
ass_expected_output$site_seq <- c('VLIPTDNsTPHKEDL', 'PVTSGLPlFFVITVT', NA)

fasta <- system.file("testdata", "reference.fasta",
                     package = "camprotR")

# Objects for testing add_peptide_positions()
app_input <- data.frame('protein'=c('L0R819', 'P98196', 'P14618'),
                        'Sequence'=c('GTRPEDSSVLIPTDNSTPHK',
                                     'QLWPTATER',
                                     'LNFSHGTHEYHAETIK'))

app_expected_output <- app_input
app_expected_output$peptide_start <- c(5, 1099, NA)
app_expected_output$peptide_end <- c(24, 1107, NA)

#### Tests ---------------------------------------------------------------------
test_that("parse_PTM_scores works", {
  expect_equal_to_reference(
    parse_PTM_scores(phospho),
    file='reference/parsed_psm_tmt_phospho.rds'
  )
})

test_that("add_PTM_positions works", {
  nucleolin_parsed <- parse_PTM_scores(nucleolin_psm)
  expect_equal_to_reference(
    add_PTM_positions(nucleolin_parsed, system.file("testdata", "nucleolin.fasta",
                                                    package = "camprotR")),
    file = "reference/parsed_nucleolin_PTM.rds"
  )
})

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
test_that("parse_PTM_scores produces error if input is incorrect class", {
  expect_error(
    parse_PTM_scores("banana"),
    "'obj' must be a data.frame"
  )
})

test_that("get_sequence produces NA for protein not in proteome", {
  expect_equal(get_sequence(proteome, 'A0R819', '1', pad=7), NA)
})

test_that("get_sequence produces NA for protein with multiple PTMs", {
  expect_equal(get_sequence(proteome, 'L0R819', '40; 120', pad=7), NA)
})

test_that("get_sequence warns if the PTM position is outside protein length", {
  expect_warning(get_sequence(proteome, 'L0R819', '200', pad=7))
})
