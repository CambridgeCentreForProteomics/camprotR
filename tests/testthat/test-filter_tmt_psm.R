context('Filter TMT PSMs')

test_that("Filter TMT PSMs", {
  data('tmt_psm', package='camprotR')
  filtered <-  filter_TMT_PSMs(tmt_psm, inter_thresh=0, sn_thresh=10)
  expect_equal_to_reference(
    MSnbase::exprs(filtered), file='filtered_tmt_psm_e.rds')
  expect_equal_to_reference(
    MSnbase::fData(filtered), file='filtered_tmt_psm_f.rds')
})
