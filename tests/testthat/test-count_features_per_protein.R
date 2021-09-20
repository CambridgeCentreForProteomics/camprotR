save_feature_count <- function(msnset, master_prot_col) {
  path <- tempfile(tmpdir = test_path("testdata"), fileext = ".txt")
  out <- count_features_per_protein(obj = msnset, master_prot_col = master_prot_col)
  write.table(out, file = path,
              sep = "\t", row.names = FALSE, col.names = TRUE)

  return(path)
}

test_that("count_features_per_protein() works", {
  load(test_path("testdata/small_psm_tmt_total.rda"))

  expect_snapshot_file(
    save_feature_count(
      msnset = small_psm_tmt_total,
      master_prot_col = "Master.Protein.Accessions"
    ),
    "feature_count.txt"
  )
})
