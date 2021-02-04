context("qc_plots")

#### Setup ---------------------------------------------------------------------
data("tmt_psm", package = "camprotR")

#### Tests ---------------------------------------------------------------------
test_that("plot_quant boxplot works", {
  expect_equal_to_reference(
    plot_quant(tmt_psm, method = "box", facet_by_sample = FALSE),
    "reference/plot_quant_box.rds"
  )
})

test_that("plot_quant histogram works", {
  expect_equal_to_reference(
    plot_quant(tmt_psm, method = "histogram", facet_by_sample = FALSE),
    "reference/plot_quant_hist.rds"
  )
})

test_that("plot_quant density works", {
  expect_equal_to_reference(
    plot_quant(tmt_psm, method = "density", facet_by_sample = FALSE),
    "reference/plot_quant_dens.rds"
  )
})

test_that("plot_quant facet_by_sample works", {
  expect_equal_to_reference(
    plot_quant(tmt_psm, method = "density", facet_by_sample = TRUE),
    "reference/plot_quant_facet.rds"
  )
})

#### Sanity checks -------------------------------------------------------------
test_that("plot_quant throws error if method is nonsense", {
  expect_error(
    plot_quant(tmt_psm, method = "banana", facet_by_sample = FALSE)
  )
})
