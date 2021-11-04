#### Setup ---------------------------------------------------------------------
library(dplyr)

#### Tests ---------------------------------------------------------------------

snap_parse_PTM_scores <- function(df) {
  # make a tempfile to compare against snapshot file
  path <- tempfile(tmpdir = test_path("testdata"), fileext = ".txt")

  # run function to test and capture output
  out <- parse_PTM_scores(
    obj = df,
    threshold = 95,
    ptm_col = "PhosphoRS.Best.Site.Probabilities",
    prob_split = "; |: ",
    collapse_delimiter = ";",
    verbose = TRUE
  )

  # save output to tempfile
  write.table(out, file = path,
              sep = "\t", row.names = FALSE, col.names = TRUE)

  # return path of tempfile
  return(path)
}

test_that("parse_PTM_scores() works", {
  # load in test MSnSet
  load(test_path("testdata/small_psm_tmt_phospho.rda"))

  # compare tempfile to snapshot file
  expect_snapshot_file(
    snap_parse_PTM_scores(
      df = MSnbase::ms2df(small_psm_tmt_phospho) # convert MSnSet to data.frame
    ),
    "psm_tmt_phospho_parsed.txt"
  )

  # test that error is produced if input is incorrect class
  expect_error(
    parse_PTM_scores(obj = "banana"),
    "'obj' must be a data.frame"
  )
})

snap_add_PTM_positions <- function(df, crap_fasta) {
  # make a tempfile to compare against snapshot file
  path <- tempfile(tmpdir = test_path("testdata"), fileext = ".txt")

  # run function to test and capture output
  out <- add_PTM_positions(
    obj = df,
    proteome_fasta = crap_fasta,
    master_protein_col = "Master.Protein.Accessions",
    sequence_col = "Sequence"
  )

  # save output to tempfile
  write.table(out, file = path,
              sep = "\t", row.names = FALSE, col.names = TRUE)

  # return path of tempfile
  return(path)
}

test_that("add_PTM_positions() works", {
  # load in test MSnSet
  load(test_path("testdata/small_psm_tmt_phospho.rda"))

  # convert MSnSet to data.frame and select a single row
  crap_ptm <- MSnbase::ms2df(small_psm_tmt_phospho) %>%
    filter(Master.Protein.Accessions == "cRAP013")

  # need to run parse_PTM_scores() first
  target <- parse_PTM_scores(crap_ptm)

  # compare tempfile to snapshot file
  expect_snapshot_file(
    snap_add_PTM_positions(
      df = target,
      crap_fasta = system.file("extdata", "cRAP_20190401.fasta.gz", package = "camprotR")
    ),
    "psm_tmt_phospho_scored.txt"
  )
})

test_that("get_sequence() works", {
  crap_proteome <- Biostrings::readAAStringSet(
    filepath = system.file("extdata", "cRAP_20190401.fasta.gz", package = "camprotR")
  )
  names(crap_proteome) <- sapply(strsplit(names(crap_proteome), split='\\|'), '[[', 3)

  # test that get_sequence() works for a protein with single PTM
  expect_equal(
    get_sequence(proteome = crap_proteome, protein = "P00330", ptm_position = 8, pad = 7),
    "MSIPETQkGVIFYES"
  )

  # test that get_sequence() can add padding at the start
  expect_equal(
    get_sequence(proteome = crap_proteome, protein = "P00330", ptm_position = 4, pad = 7),
    "___MSIPeTQKGVI"
  )

  # test that get_sequence() can add padding at the end
  expect_equal(
    get_sequence(proteome = crap_proteome, protein = "P00330", ptm_position = 346, pad = 7),
    "VGRYVVDtSK_____"
  )

  # test that get_sequence() produces NA for protein not in proteome
  expect_equal(
    get_sequence(proteome = crap_proteome, protein = "AAAAAA", ptm_position = 346, pad = 7),
    NA
  )

  # test that get_sequence() produces NA for protein with multiple PTMs
  expect_equal(
    get_sequence(proteome = crap_proteome, protein = "P00330", ptm_position = "40; 120", pad = 7),
    NA
  )

  # test that get_sequence() warns if PTM position is outside protein length
  expect_warning(
    get_sequence(proteome = crap_proteome, protein = "P00330", ptm_position = 1000, pad = 7)
  )
})

test_that("add_site_sequence() works", {
  target <- data.frame(
    Master.Protein.Accessions = c("cRAP001", "cRAP002", "XXXXXX"),
    ptm_position = c(10, 28, 40) # column must have this name
  )

  expect_equal(
    add_site_sequence(
      obj = target,
      proteome_fasta = system.file("extdata", "cRAP_20190401.fasta.gz", package = "camprotR"),
      master_protein_col = "Master.Protein.Accessions"
    )$site_seq,
    c("IPETQKGvIFYESHG", "VFRRDAHkSEVAHRF", NA)
  )
})

test_that("add_peptide_positions() works", {
  target <- data.frame(
    Master.Protein.Accessions = c("cRAP001", "cRAP002", "XXXXXX"),
    Sequence = c("LAQVAPILCAGITVYKA",
                 "FVEKCCKADDKETCFA",
                 "ACDEPQRWY")
  )

  # find peptide positions in FASTA and add peptide_start and peptide_end columns
  target <- add_peptide_positions(
    obj = target,
    proteome_fasta = system.file("extdata", "cRAP_20190401.fasta.gz", package = "camprotR"),
    master_protein_col = "Master.Protein.Accessions"
  )

  # check peptide_start is of expected values
  expect_equal(
    target$peptide_start,
    c(146, 578, NA)
  )

  # check peptide_end is of expected values
  expect_equal(
    target$peptide_end,
    c(162, 593, NA)
  )
})
