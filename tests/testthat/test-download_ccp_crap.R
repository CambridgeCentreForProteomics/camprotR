test_that("check_uniprot_release() works", {
  expect_true(grepl("[0-9]{4}_[0-9]{2}", check_uniprot_release()))
})

test_that("sub_crap() works", {
  expect_equal(sub_crap("|text|"), "|cRAP001|text|")
  expect_equal(sub_crap("|text|", 77), "|cRAP077|text|")
  expect_equal(sub_crap("|text|", 111, 4), "|cRAP0111|text|")
})

test_that("download_ccp_crap() works", {
  withr::with_tempfile("tf", {
    download_ccp_crap(
      file = tf,
      is_crap = TRUE,
      overwrite = FALSE,
      verbose = TRUE
    )

    target <- Biostrings::readAAStringSet(filepath = tf)

    # Use expect_snapshot_output(str()) rather than expect_snapshot_file(), which
    # should hopefully be more robust to small changes in UniProt headers or
    # sequences that can occur over time.
    expect_snapshot_output(str(target@ranges))
  })
})
