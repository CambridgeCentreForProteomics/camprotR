context("tmt_qc_plots")

#### Setup ---------------------------------------------------------------------

#### Tests ---------------------------------------------------------------------
test_that("Tallying sub-notch values in PSM-level TMT data, extract info", {
  data('psm_data_notch_count', package='camprotR')
  expect_equal_to_reference(get_notch_per_protein(psm_data_notch_count), 'reference/notch_count.rds')
})

#### Sanity checks -------------------------------------------------------------
