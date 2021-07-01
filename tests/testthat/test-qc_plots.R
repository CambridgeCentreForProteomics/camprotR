context("qc_plots")

#### Setup ---------------------------------------------------------------------
load(system.file("testdata", "small_psm_tmt_total.rda",
                 package = "camprotR"))

#### Tests ---------------------------------------------------------------------
# test_that("plot_quant boxplot works", {
#   expect_equal_to_reference(
#     plot_quant(small_psm_tmt_total, method = "box", facet_by_sample = FALSE),
#     "reference/plot_quant_box.rds"
#   )
# })
#
# test_that("plot_quant histogram works", {
#   expect_equal_to_reference(
#     plot_quant(small_psm_tmt_total, method = "histogram", facet_by_sample = FALSE),
#     "reference/plot_quant_hist.rds"
#   )
# })
#
# test_that("plot_quant density works", {
#   expect_equal_to_reference(
#     plot_quant(small_psm_tmt_total, method = "density", facet_by_sample = FALSE),
#     "reference/plot_quant_dens.rds"
#   )
# })
#
# test_that("plot_quant facet_by_sample works", {
#   expect_equal_to_reference(
#     plot_quant(small_psm_tmt_total, method = "density", facet_by_sample = TRUE),
#     "reference/plot_quant_facet.rds"
#   )
# })

#### Sanity checks -------------------------------------------------------------
test_that("plot_quant throws error if method is nonsense", {
  expect_error(
    plot_quant(small_psm_tmt_total, method = "banana", facet_by_sample = FALSE)
  )
})
