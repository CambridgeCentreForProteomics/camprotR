#### Functions to create test objects ------------------------------------------

create_small_psm_tmt_phospho <- function(df) {
  # set seed for reproducible sampling
  set.seed(2021)

  # randomly sample 234 rows
  df_filt <- df[sample(nrow(df), 234), ]

  # extract expression data (i.e. PSM abundances in each TMT channel)
  abundance_cols <- colnames(df_filt)[grepl('Abundance.', colnames(df_filt))]
  e_data <- as.matrix(df_filt[, abundance_cols])

  # extract feature data that we want (i.e. metadata about each PSM)
  f_data <- df_filt[, c("Master.Protein.Accessions"), drop = FALSE]

  # construct MSnSet
  small_psm_tmt_phospho <- MSnbase::MSnSet(exprs = e_data, fData = f_data)

  # output .rda file
  save(small_psm_tmt_phospho, file = test_path("testdata/small_psm_tmt_phospho.rda"))
}

create_small_psm_tmt_total <- function(df) {
  # filter for first 100 PSMs
  df_filt <- df[1:100, ]

  # extract expression data (i.e. PSM abundances in each TMT channel)
  abundance_cols <- colnames(df_filt)[grepl('Abundance.', colnames(df_filt))]
  e_data <- as.matrix(df_filt[, abundance_cols])

  # extract feature data that we want (i.e. metadata about each PSM)
  f_data <- df_filt[, c("Master.Protein.Accessions",
                        "Isolation.Interference.in.Percent",
                        "Average.Reporter.SN")]

  # construct MSnSet
  small_psm_tmt_total <- MSnbase::MSnSet(exprs = e_data, fData = f_data)

  # output .rda file
  save(small_psm_tmt_total, file = test_path("testdata/small_psm_tmt_total.rda"))
}

create_small_pep_silac_p0 <- function(df) {
  # select only necessary columns
  small_pep_silac_p0 <- df %>%
    select(
      Sequence,
      Modifications,
      Master.Protein.Accessions,
      Protein.Accessions,
      Number.of.Protein.Groups,
      Quan.Info,
      Abundances.Grouped.F1.Light,
      Abundances.Grouped.F1.Heavy
    )

  # output .txt file
  write.table(small_pep_silac_p0, file = test_path("testdata/small_pep_silac_p0.txt"),
              sep = "\t", row.names = FALSE, col.names = TRUE)
}

create_small_pep_silac_p4 <- function(df) {
  # select only necessary columns
  small_pep_silac_p4 <- df %>%
    select(
      Sequence,
      Modifications,
      Master.Protein.Accessions,
      Protein.Accessions,
      Number.of.Protein.Groups,
      Quan.Info,
      Abundances.Grouped.F5.Light,
      Abundances.Grouped.F5.Heavy
    )

  # output .txt file
  write.table(small_pep_silac_p4, file = test_path("testdata/small_pep_silac_p4.txt"),
              sep = "\t", row.names = FALSE, col.names = TRUE)
}

#### Functions to setup and cleanup the test environment -----------------------

setup_testenv <- function() {
  # make temp folder to store data used across multiple tests
  dir.create("testdata")

  # create a random file to test setup has worked (see test-setup.R)
  write.table(data.frame(letters = c("A", "B", "C")),
              file = test_path("testdata/letters.txt"),
              sep = "\t", row.names = FALSE, col.names = TRUE)

  # create small MSnSet from psm_tmt_phospho for testing purposes
  create_small_psm_tmt_phospho(df = psm_tmt_phospho)

  # create small MSnSet from psm_tmt_total for testing purposes
  create_small_psm_tmt_total(df = psm_tmt_total)

  # create small .txt file from pep_silac_p0 for testing purposes
  create_small_pep_silac_p0(df = pep_silac_p0)

  # create small .txt file from pep_silac_p4 for testing purposes
  create_small_pep_silac_p4(df = pep_silac_p4)
}

cleanup_testenv <- function() {
  # remove testdata temp folder and contents
  unlink("testdata", recursive = TRUE)
}

#### Actually setup and cleanup the test environment ---------------------------

# setup test environment
setup_testenv()

# ... tests happen here ...

# cleanup test environment when finished
withr::defer(cleanup_testenv(), testthat::teardown_env())
