context('normalise')

#### Setup ---------------------------------------------------------------------
load(system.file("tests", "testthat", "data-test", "small_psm_tmt_total.rda",
                 package = "camprotR"))

#### Tests ---------------------------------------------------------------------
test_that("Center-normalise, no reference", {
  expect_equal_to_reference(center_normalise_to_ref(small_psm_tmt_total,
                                                    get_medians(small_psm_tmt_total)),
                            file='reference/tmt_center_normalised.rds')
})

#### Sanity checks -------------------------------------------------------------
