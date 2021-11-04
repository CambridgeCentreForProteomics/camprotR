snap_count_features <- function(msnset, master_prot_col) {
  # make a tempfile to compare against snapshot file
  path <- tempfile(tmpdir = test_path("testdata"), fileext = ".txt")

  # run function to test and capture output
  out <- count_features_per_protein(obj = msnset, master_prot_col = master_prot_col)

  # save output to tempfile
  write.table(out, file = path,
              sep = "\t", row.names = FALSE, col.names = TRUE)

  # return path of tempfile
  return(path)
}

test_that("count_features_per_protein() works", {
  # load in test MSnSet
  load(test_path("testdata/small_psm_tmt_total.rda"))

  # compare tempfile to snapshot file
  expect_snapshot_file(
    snap_count_features(
      msnset = small_psm_tmt_total,
      master_prot_col = "Master.Protein.Accessions"
    ),
    "feature_count.txt"
  )
})
