context('count_features_per_protein')

#### Setup ---------------------------------------------------------------------
load(system.file("testdata", "small_psm_tmt_total.rda",
                 package = "camprotR"))

#### Tests ---------------------------------------------------------------------
test_that("Count PSMs", {
  expect_equal_to_reference(count_features_per_protein(small_psm_tmt_total),
                            file='reference/tmt_feature_counts.rds')
})

#### Sanity checks -------------------------------------------------------------
