test_that("append_fasta() works", {
  # create self cleaning tempfile
  withr::with_tempfile("tf", {
    # make some fake fasta
    Biostrings::writeXStringSet(
      x = Biostrings::AAStringSet(x = c("a fasta header" = "QQQPPPMMM")),
      filepath = tf
    )

    # append 'commercial_reagents' fasta to fake fasta file
    append_fasta(
      file1 = system.file("extdata/commercial_reagents.fasta", package = "camprotR"),
      file2 = tf,
      is_crap = TRUE,
      crap_start = 1
    )

    # read in tempfile which should have been modified by append_fasta()
    target <- Biostrings::readAAStringSet(tf)

    # check fasta lengths are correct
    expect_equal(lengths(target), c(9, 274, 261))

    # check fasta headers are okay
    expect_true(grepl("Endoproteinase GluC", names(target)[2]))
    expect_true(grepl("rLys-C", names(target)[3]))

    # check CRAP number has been added
    expect_true(grepl("cRAP001", names(target)[2]))
    expect_true(grepl("cRAP002", names(target)[3]))
  })
})
