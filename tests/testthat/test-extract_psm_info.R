context('extract_psm_info')

#### Setup ---------------------------------------------------------------------

load(system.file("testdata", "small_psm_silac.rda", package = "camprotR"))

#### Tests ---------------------------------------------------------------------
test_that("silac_psm_seq_int works", {
  expect_equal_to_reference(
    silac_psm_seq_int(small_psm_silac, sequence_col='Annotated.Sequence'),
    "reference/silac_psm_seq_int.rds"
  )
})

#### Sanity checks -------------------------------------------------------------
