snap_center_normalise_to_ref <- function(msnset) {
  # make a tempfile to compare against snapshot file
  path <- tempfile(tmpdir = test_path("testdata"), fileext = ".txt")

  # run function to test and capture output
  out <- center_normalise_to_ref(
    obj = msnset,
    medians = get_medians(obj = msnset)
  )

  # save output to tempfile
  write.table(exprs(out), file = path,
              sep = "\t", row.names = FALSE, col.names = TRUE)

  # return path of tempfile
  return(path)
}

test_that("center_normalise_to_ref() works", {
  # load in test MSnSet
  load(test_path("testdata/small_psm_tmt_total.rda"))

  # compare tempfile to snapshot file
  expect_snapshot_file(
    snap_center_normalise_to_ref(
      msnset = small_psm_tmt_total
    ),
    "center_normalise_to_ref.txt",
    cran = TRUE
  )
})
