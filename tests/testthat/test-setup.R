test_that("self cleaning test infrastructure works", {
  df <- read.delim(file = test_path("testdata/letters.txt"), header = TRUE)
  expect_equal(df$letters, c("A", "B", "C"))
})
