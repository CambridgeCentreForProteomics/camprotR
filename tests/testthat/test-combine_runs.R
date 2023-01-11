snap_get_parsimony_pep2prot <- function(infiles) {
  # make a tempfile to compare against snapshot file
  path <- tempfile(tmpdir = test_path("testdata"), fileext = ".txt")

  # run function to test and capture output
  out <- get_parsimony_pep2prot(
    infiles = infiles,
    seq_col = "Sequence",
    prot_col = "Protein.Accessions",
    master_prot_col = "Master.Protein.Accessions",
    compare_old_new = TRUE
  )

  # save output to tempfile
  write.table(out, file = path,
              sep = "\t", row.names = FALSE, col.names = TRUE)

  # return path of tempfile
  return(path)
}

test_that("get_parsimony_pep2prot() works", {
  # compare tempfile to snapshot file
  expect_snapshot_file(
    snap_get_parsimony_pep2prot(
      list(
        test_path("testdata/small_pep_silac_p0.txt"),
        test_path("testdata/small_pep_silac_p4.txt")
      )
    ),
    "parsimony.txt",
    cran = TRUE
  )
})
