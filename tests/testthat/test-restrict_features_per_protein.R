context("restrict_features_per_protein")

#### Setup ---------------------------------------------------------------------
# Load small PSM TMT MSnSet
load(system.file("testdata", "small_psm_tmt_total.rda",
                 package = "camprotR"))

# Filter the PSMs to only retain those with S:N > 10 using `filter_TMT_PSMs`
# Using the same function, we also remove PSMs with interference/co-isolation > 50%
psm_filt <- filter_TMT_PSMs(small_psm_tmt_total, inter_thresh = 50, sn_thresh = 10)

#### Tests ---------------------------------------------------------------------
test_that("restrict_features_per_protein works", {
  expect_equal_to_reference(
    MSnbase::exprs(restrict_features_per_protein(psm_filt, min_features = 2)),
    "reference/restrict_features.rds"
  )
})

#### Sanity checks -------------------------------------------------------------
