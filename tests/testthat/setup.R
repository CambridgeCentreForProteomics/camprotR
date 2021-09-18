setup_testenv <- function() {
  # make temp folder to store data used across multiple tests
  dir.create("testdata")

  # create a random file to test setup has worked (see test-setup.R)
  write.table(data.frame(letters = c("A", "B", "C")), file = test_path("testdata/letters.txt"),
              sep = "\t", row.names = FALSE, col.names = TRUE)
}

cleanup_testenv <- function() {
  # remove testdata temp folder and contents
  unlink("testdata", recursive = TRUE)
}

# setup test environment
setup_testenv()

# ... tests happen here ...

# cleanup test environment when finished
withr::defer(cleanup_testenv(), testthat::teardown_env())
