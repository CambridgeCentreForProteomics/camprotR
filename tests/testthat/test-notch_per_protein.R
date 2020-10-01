context('Tallying sub-notch values in PSM-level TMT data')

test_that("Extract info", {
  data('psm_data_notch_count', package='camprotR')
  expect_equal_to_reference(get_notch_per_protein(psm_data_notch_count), 'notch_count.rds')
})


