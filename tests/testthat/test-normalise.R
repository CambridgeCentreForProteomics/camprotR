save_center_normalise_to_ref <- function(msnset) {
  path <- tempfile(tmpdir = test_path("testdata"), fileext = ".txt")
  out <- center_normalise_to_ref(
    obj = msnset,
    medians = get_medians(obj = msnset)
  )
  write.table(exprs(out), file = path,
              sep = "\t", row.names = FALSE, col.names = TRUE)

  return(path)
}

test_that("center_normalise_to_ref() works", {
  load(test_path("testdata/small_psm_tmt_total.rda"))

  expect_snapshot_file(
    save_center_normalise_to_ref(
      msnset = small_psm_tmt_total
    ),
    "center_normalise_to_ref.txt"
  )
})
