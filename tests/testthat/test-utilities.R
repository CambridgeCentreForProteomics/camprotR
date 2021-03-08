context("utilities")

#### Setup ---------------------------------------------------------------------

#### Tests ---------------------------------------------------------------------

test_that("remove_dots works", {
  df <- data.frame(
    column...name = c(1, 2, 3)
  )

  colnames(df) <- remove_dots(colnames(df))
  expect_equal(grep("\\.", colnames(df)), 1)
})
