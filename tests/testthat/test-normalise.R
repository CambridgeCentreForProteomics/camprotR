context('Center-normalise')

test_that("Center-normalise, no reference", {
  data('tmt_psm', package='camprotR')
  expect_equal_to_reference(center_normalise_to_ref(tmt_psm, get_medians(tmt_psm)),
                            file='tmt_center_normalised.rds')
})
