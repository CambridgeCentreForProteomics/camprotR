test_that("plot_quant() works", {
  # load in test MSnSet
  load(test_path("testdata/small_psm_tmt_total.rda"))

  # filter out NA
  target <- MSnbase::filterNA(small_psm_tmt_total)

  # boxplot
  vdiffr::expect_doppelganger(
    "plot-quant-method-box",
    plot_quant(target, method = 'box', facet_by_sample = FALSE)
  )

  # histogram
  vdiffr::expect_doppelganger(
    "plot-quant-method-histogram",
    plot_quant(target, method = 'histogram', facet_by_sample = FALSE)
  )

  # density
  vdiffr::expect_doppelganger(
    "plot-quant-method-density",
    plot_quant(target, method = 'density', facet_by_sample = FALSE)
  )

  # facet_by_sample
  vdiffr::expect_doppelganger(
    "plot-quant-method-density-facet",
    plot_quant(target, method = 'density', facet_by_sample = TRUE)
  )

  # error occurs if method is nonsense
  expect_error(
    plot_quant(target, method = 'banana', facet_by_sample = FALSE)
  )
})
