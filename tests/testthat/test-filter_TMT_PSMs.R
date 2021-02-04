context('filter_TMT_PSMs')

#### Setup ---------------------------------------------------------------------

#### Tests ---------------------------------------------------------------------
test_that("Filter TMT PSMs", {
  data('tmt_psm', package='camprotR')
  filtered <-  filter_TMT_PSMs(tmt_psm, inter_thresh=0, sn_thresh=10)
  expect_equal_to_reference(
    MSnbase::exprs(filtered), file='reference/filtered_tmt_psm_e.rds')
  expect_equal_to_reference(
    MSnbase::fData(filtered), file='reference/filtered_tmt_psm_f.rds')
})

test_that("Update average S/N", {
  data('tmt_psm', package='camprotR')
  filtered <-  update_average_sn(tmt_psm)
  expect_equal_to_reference(
    MSnbase::fData(filtered), file='reference/tmt_psm_f_update_sn.rds')
})

#### Sanity checks -------------------------------------------------------------
