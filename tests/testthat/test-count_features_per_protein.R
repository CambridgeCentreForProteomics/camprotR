context('count_features_per_protein')

#### Setup ---------------------------------------------------------------------

#### Tests ---------------------------------------------------------------------
test_that("Count PSMs", {
  data('tmt_psm', package='camprotR')
  expect_equal_to_reference(count_features_per_protein(tmt_psm),
                            file='reference/tmt_feature_counts.rds')
})

#### Sanity checks -------------------------------------------------------------
