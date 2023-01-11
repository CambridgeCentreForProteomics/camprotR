test_that("plot_TMT_notch() works", {
  # load in test MSnSet
  load(test_path("testdata/small_psm_tmt_phospho.rda"))

  # test that defaults work
  vdiffr::expect_doppelganger(
    "plot-tmt-notch-default",
    plot_TMT_notch(MSnbase::filterNA(small_psm_tmt_phospho))
  )

  # test that you can facet by sample
  vdiffr::expect_doppelganger(
    "plot-tmt-notch-facet",
    plot_TMT_notch(MSnbase::filterNA(small_psm_tmt_phospho), facet_by_sample = TRUE)
  )
})

snap_get_notch_per_protein <- function(msnset) {
  # make a tempfile to compare against snapshot file
  path <- tempfile(tmpdir = test_path("testdata"), fileext = ".txt")

  # run function to test and capture output
  out <- get_notch_per_protein(
    obj = msnset,
    master_prot_col = "Master.Protein.Accessions",
    notch_upper = 5.75
  )

  # save output to tempfile
  write.table(out, file = path,
              sep = "\t", row.names = FALSE, col.names = TRUE)

  # return path of tempfile
  return(path)
}

test_that("get_notch_per_protein() works", {
  # load in test MSnSet
  load(test_path("testdata/small_psm_tmt_phospho.rda"))

  # compare tempfile to snapshot file
  expect_snapshot_file(
    snap_get_notch_per_protein(
      msnset = small_psm_tmt_phospho
    ),
    "notch_per_protein.txt",
    cran = TRUE
  )
})

test_that("plot_below_notch_per_prot() works", {
  # load in test MSnSet
  load(test_path("testdata/small_psm_tmt_phospho.rda"))

  # test that defaults work
  vdiffr::expect_doppelganger(
    "plot-below-notch-default",
    plot_below_notch_per_prot(
      get_notch_per_protein(small_psm_tmt_phospho)
    )
  )
})

test_that("plot_fraction_below_notch_per_prot() works", {
  # load in test MSnSet
  load(test_path("testdata/small_psm_tmt_phospho.rda"))

  # test that defaults work
  vdiffr::expect_doppelganger(
    "plot-fraction-below-notch-default",
    plot_fraction_below_notch_per_prot(
      get_notch_per_protein(small_psm_tmt_phospho)
    )
  )
})

test_that("plot_missing_SN() works", {
  # load in test MSnSet
  load(test_path("testdata/small_psm_tmt_total.rda"))

  # test that defaults work
  vdiffr::expect_doppelganger(
    "plot-missing-sn",
    plot_missing_SN(small_psm_tmt_total)
  )
})

test_that("plot_missing_SN_per_sample() works", {
  # load in test MSnSet
  load(test_path("testdata/small_psm_tmt_total.rda"))

  # test that defaults work
  vdiffr::expect_doppelganger(
    "plot-missing-sn-per-sample",
    plot_missing_SN_per_sample(small_psm_tmt_total)
  )
})
