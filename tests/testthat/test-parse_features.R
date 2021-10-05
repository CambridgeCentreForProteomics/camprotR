#### Setup ---------------------------------------------------------------------

# get headers from CCP cRAP file
crap_fasta <- Biostrings::fasta.index(
  system.file("extdata", "cRAP_20190401.fasta.gz", package = "camprotR"),
  seqtype = "AA"
)

# extract the non cRAP UniProt accessions associated with each cRAP protein
# which we will use for filtering features in parse_features()
crap_accessions <- unlist(
  regmatches(
    crap_fasta$desc,
    gregexpr("(?<=\\|).*?(?=\\|)", crap_fasta$desc, perl = TRUE)
  )
)

#### Tests ---------------------------------------------------------------------

snap_parse_features <- function(df, silac, TMT, level, crap_proteins) {
  # make a tempfile to compare against snapshot file
  path <- tempfile(tmpdir = test_path("testdata"), fileext = ".txt")

  # run function to test and capture output
  out <- parse_features(
    data = df,
    master_protein_col = "Master.Protein.Accessions",
    protein_col = "Protein.Accessions",
    unique_master = TRUE,
    silac = silac,
    TMT = TMT,
    level = level,
    filter_crap = TRUE,
    crap_proteins = crap_proteins,
    filter_associated_crap = TRUE
  )

  # save output to tempfile
  write.table(out, file = path,
              sep = "\t", row.names = FALSE, col.names = TRUE)

  # return path of tempfile
  return(path)
}

test_that("parse_features() works with TMT PSM data", {
  # compare tempfile to snapshot file
  expect_snapshot_file(
    snap_parse_features(
      df = psm_tmt_total,
      silac = FALSE,
      TMT = TRUE,
      level = "PSM",
      crap_proteins = crap_accessions
    ),
    "psm_tmt_total_parsed.txt"
  )
})

test_that("parse_features() works with SILAC peptide data", {
  # load in test data.frame
  small_pep_silac_p4 <- read.delim(test_path("testdata/small_pep_silac_p4.txt"))

  # compare tempfile to snapshot file
  expect_snapshot_file(
    snap_parse_features(
      df = psm_tmt_total,
      silac = TRUE,
      TMT = FALSE,
      level = "peptide",
      crap_proteins = crap_accessions
    ),
    "small_pep_silac_p4_parsed.txt"
  )
})

#### Sanity checks -------------------------------------------------------------

test_that("An error is produced if crap_proteins is not given", {
  expect_error(
    parse_features(psm_data, TMT = TRUE, level = "PSM"),
    "must supply the crap_proteins argument to filter cRAP proteins"
  )
})
