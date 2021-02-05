context("tmt_qc_plots")

#### Setup ---------------------------------------------------------------------

# Load some TMT PSM level data for tests
data('psm_data_notch_count', package = 'camprotR')
data("tmt_psm", package = "camprotR")

#### Tests ---------------------------------------------------------------------
# bit of a hacky test because can't use expect_equal_to_reference as the
# object changes slightly every time you run the test.
test_that("plot_TMT_notch works", {
  p <- plot_TMT_notch(MSnbase::filterNA(psm_data_notch_count))
  expect_silent(ggplotGrob(p))
})

test_that("plot_TMT_notch faceted by sample works", {
  p <- plot_TMT_notch(MSnbase::filterNA(psm_data_notch_count), facet_by_sample = TRUE)
  expect_silent(ggplotGrob(p))
})

test_that("get_notch_per_protein works", {
  expect_equal_to_reference(
    get_notch_per_protein(psm_data_notch_count),
    'reference/notch_count.rds'
  )
})

test_that("plot_below_notch_per_prot works", {
  expect_equal_to_reference(
    plot_below_notch_per_prot(
      get_notch_per_protein(psm_data_notch_count)
    ),
    "reference/plot_below_notch_per_prot.rds"
  )
})

test_that("plot_fraction_below_notch_per_prot works", {
  expect_equal_to_reference(
    plot_fraction_below_notch_per_prot(
      get_notch_per_protein(psm_data_notch_count)
    ),
    "reference/plot_fraction_below_notch_per_prot.rds"
  )
})

test_that("plot_missing_SN works", {
  expect_equal_to_reference(
    plot_missing_SN(tmt_psm),
    "reference/plot_missing_SN.rds"
  )
})

test_that("plot_missing_SN_per_sample works", {
  expect_equal_to_reference(
    plot_missing_SN_per_sample(tmt_psm),
    "reference/plot_missing_SN_per_sample.rds"
  )
})

#### Sanity checks -------------------------------------------------------------
