test_that("self cleaning test infrastructure works", {
  df <- read.delim(file = test_path("testdata/letters.txt"), header = TRUE)
  expect_equal(df$letters, c("A", "B", "C"))

  load(test_path("testdata/small_psm_tmt_phospho.rda"))
  expect_equal(dim(small_psm_tmt_phospho), c(234, 10))

  load(test_path("testdata/small_psm_tmt_total.rda"))
  expect_equal(dim(small_psm_tmt_total), c(100, 10))

  small_pep_silac_p0 <- read.delim(test_path("testdata/small_pep_silac_p0.txt"))
  expect_equal(dim(small_pep_silac_p0), c(418, 8))

  small_pep_silac_p4 <- read.delim(test_path("testdata/small_pep_silac_p4.txt"))
  expect_equal(dim(small_pep_silac_p4), c(470, 8))
})
