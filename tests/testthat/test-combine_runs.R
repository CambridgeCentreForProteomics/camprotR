save_parsimony <- function(infiles) {
  path <- tempfile(tmpdir = test_path("testdata"), fileext = ".txt")
  out <- get_parsimony_pep2prot(
    infiles = infiles,
    seq_col = "Sequence",
    prot_col = "Protein.Accessions",
    master_prot_col = "Master.Protein.Accessions",
    compare_old_new = TRUE
  )
  write.table(out, file = path,
              sep = "\t", row.names = FALSE, col.names = TRUE)

  return(path)
}

test_that("get_parsimony_pep2prot() works", {
  expect_snapshot_file(
    save_parsimony(
      list(
        test_path("testdata/small_pep_silac_p0.txt"),
        test_path("testdata/small_pep_silac_p4.txt")
      )
    ),
    "parsimony.txt"
  )
})
