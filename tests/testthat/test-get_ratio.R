context("get_ratio")

#### Setup ---------------------------------------------------------------------

df <- data.frame(
  treatment = c(30, NA, 12, NA),
  control = c(20, 12, NA, NA)
)

#### Tests ---------------------------------------------------------------------

test_that("get_ratio outputs a data.frame by default", {
  df2 <- get_ratio(df, treatment, control)

  expect_s3_class(df2, "data.frame")
})

test_that("get_ratio outputs a named list if bind = FALSE", {
  df2 <- get_ratio(df, treatment, control, bind = FALSE)

  expect_type(df2, "list")
  expect_named(df2)
})

test_that("there are no NAs in the new 'missing' column", {
  df2 <- get_ratio(df, treatment, control)

  expect_false(any(is.na(df2$missing)))
})
