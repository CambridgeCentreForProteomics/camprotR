context("test_function")

#### Setup ---------------------------------------------------------------------

df <- data.frame(
  x = c(0, 1, 2),
  y = c(0, 1, 2)
)

#### Tests ---------------------------------------------------------------------

test_that("test_function can create a ggplot", {
  plot1 <- test_function(df, x, y)

  expect_length(layer_grob(plot1), 1)
})
