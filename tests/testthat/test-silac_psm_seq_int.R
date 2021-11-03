snap_psm_seq_int <- function(infile, sequence_col) {
  # make a tempfile to compare against snapshot file
  path <- tempfile(tmpdir = test_path("testdata"), fileext = ".txt")

  # read in data to for function to use
  df <- read.table(file = infile, header = TRUE, sep = "\t")

  # run function to test and capture output
  out <- silac_psm_seq_int(df, sequence_col = sequence_col)

  # save output to tempfile
  write.table(out, file = path,
              sep = "\t", row.names = FALSE, col.names = TRUE)

  # return path of tempfile
  return(path)
}

test_that("silac_psm_seq_int() works", {
  # compare tempfile to snapshot file
  expect_snapshot_file(
    snap_psm_seq_int(
      infile = test_path("testdata/small_psm_silac_p4.txt"),
      sequence_col = "Annotated.Sequence"
    ),
    "silac_psm_seq_int.txt"
  )
})
