test_that("append_fasta() works", {
  withr::with_tempfile(c("tf1", "tf2"), {
    # save FUBP3 fasta to tempfile 1
    httr::GET("https://www.uniprot.org/uniprot/Q96I24.fasta", config = httr::write_disk(path = tf1))

    # save GAPDH fasta to tempfile 2
    httr::GET("https://www.uniprot.org/uniprot/P04406.fasta", config = httr::write_disk(path = tf2))

    # append FUBP3 fasta to GAPDH fasta file
    append_fasta(
      file1 = tf1,
      file2 = tf2,
      is_crap = TRUE,
      crap_start = 1
    )

    # read in tempfile2 which should have been modified by append_fasta
    target <- Biostrings::readAAStringSet(tf2)

    # check fasta lengths are correct
    expect_equal(lengths(target), c(335, 572))

    # check fasta headers are okay
    expect_true(grepl("P04406", names(target)[1]))
    expect_true(grepl("Q96I24", names(target)[2]))

    # check CRAP number has been added
    expect_true(grepl("cRAP001", names(target)[2]))
  })
})
