save_psm_seq_int <- function(infile, sequence_col) {
  path <- tempfile(tmpdir = test_path("testdata"), fileext = ".txt")

  df <- read.table(file = infile, header = TRUE, sep = "\t")
  out <- silac_psm_seq_int(df, sequence_col = sequence_col)
  write.table(out, file = path,
              sep = "\t", row.names = FALSE, col.names = TRUE)

  return(path)
}

test_that("silac_psm_seq_int() works", {

  expect_snapshot_file(
    save_psm_seq_int(
      infile = test_path("testdata/small_psm_silac_p4.txt"),
      sequence_col = "Annotated.Sequence"
    ),
    "silac_psm_seq_int.txt"
  )
})
