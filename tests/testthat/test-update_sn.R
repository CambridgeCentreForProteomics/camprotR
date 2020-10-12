context('Update average S/N')

test_that("Filter TMT PSMs", {
  data('tmt_psm', package='camprotR')
  filtered <-  update_average_sn(tmt_psm)
  expect_equal_to_reference(
    MSnbase::fData(filtered), file='tmt_psm_f_update_sn.rds')
})
