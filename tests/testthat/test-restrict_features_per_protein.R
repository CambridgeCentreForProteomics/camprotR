snap_restrict_features <- function(msnset, min_features, master_protein_col, plot) {
  # make a tempfile to compare against snapshot file
  path <- tempfile(tmpdir = test_path("testdata"), fileext = ".txt")

  # run function to test and capture output
  out <- restrict_features_per_protein(
    obj = msnset,
    min_features = min_features,
    master_protein_col = master_protein_col,
    plot = plot
  )

  if (plot) {
    return(out)
  } else {
    # save output to tempfile
    write.table(MSnbase::ms2df(out), file = path,
                sep = "\t", row.names = FALSE, col.names = TRUE)

    # return path of tempfile
    return(path)
  }
}

test_that("restrict_features_per_protein() works", {
  # load in test MSnSet
  load(test_path("testdata/small_psm_tmt_total.rda"))

  target <- filter_TMT_PSMs(small_psm_tmt_total, inter_thresh = 50, sn_thresh = 10)

  # test that plot = TRUE works
  vdiffr::expect_doppelganger(
    "restrict-features-histogram",
    snap_restrict_features(
      msnset = target,
      min_features = 2,
      master_protein_col = "Master.Protein.Accessions",
      plot = TRUE
    )
  )

  # test that plot = FALSE works
  expect_snapshot_file(
    snap_restrict_features(
      msnset = target,
      min_features = 2,
      master_protein_col = "Master.Protein.Accessions",
      plot = FALSE
    ),
    "restrict-features-msnset.txt",
    cran = TRUE
  )
})
