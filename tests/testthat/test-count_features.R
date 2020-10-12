context('Count features per protein')

test_that("Filter TMT PSMs", {
  data('tmt_psm', package='camprotR')
  expect_equal_to_reference(count_features_per_protein(tmt_psm),
                            file='tmt_feature_counts.rds')
})
