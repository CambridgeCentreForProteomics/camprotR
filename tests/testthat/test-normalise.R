context('normalise')

#### Setup ---------------------------------------------------------------------

#### Tests ---------------------------------------------------------------------
test_that("Center-normalise, no reference", {
  data('tmt_psm', package='camprotR')
  expect_equal_to_reference(center_normalise_to_ref(tmt_psm, get_medians(tmt_psm)),
                            file='reference/tmt_center_normalised.rds')
})

#### Sanity checks -------------------------------------------------------------
