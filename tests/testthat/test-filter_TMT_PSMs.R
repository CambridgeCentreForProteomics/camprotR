context('filter_TMT_PSMs')

#### Setup ---------------------------------------------------------------------
load(system.file("tests", "testthat", "data-test", "small_psm_tmt_total.rda",
                 package = "camprotR"))

#### Tests ---------------------------------------------------------------------
test_that("Filter TMT PSMs", {
  filtered <-  filter_TMT_PSMs(small_psm_tmt_total, inter_thresh=0, sn_thresh=10)
  expect_equal_to_reference(
    MSnbase::exprs(filtered), file='reference/filtered_tmt_psm_e.rds')
  expect_equal_to_reference(
    MSnbase::fData(filtered), file='reference/filtered_tmt_psm_f.rds')
})

test_that("Update average S/N", {
  filtered <-  update_average_sn(small_psm_tmt_total)
  expect_equal_to_reference(
    MSnbase::fData(filtered), file='reference/tmt_psm_f_update_sn.rds')
})

#### Sanity checks -------------------------------------------------------------
