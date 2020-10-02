context('camprotR theme')
# setup #
p <- ggplot(data.frame('x'=1:10, 'y'=1:10), aes(x, y)) + geom_point()
#########

test_that('default theme', code={
  expect_equal_to_reference(theme_camprot(), 'camprotR_theme.rds')
})

test_that('base_size=10', code={
  expect_equal_to_reference(theme_camprot(base_size=10), 'camprotR_theme_10.rds')
})

test_that('default theme', code={
  expect_equal_to_reference(theme_camprot(base_family='sans'), 'camprotR_theme_sans.rds')
})

test_that('default theme', code={
  expect_equal_to_reference(theme_camprot(aspect_square=FALSE), 'camprotR_theme_not_square.rds')
})

test_that('default theme', code={
  expect_equal_to_reference(theme_camprot(border=FALSE), 'camprotR_theme_no_border.rds')
})
