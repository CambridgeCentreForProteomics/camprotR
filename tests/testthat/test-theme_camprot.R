context('theme_camprot')

#### Setup ---------------------------------------------------------------------
p <- ggplot(data.frame('x'=1:10, 'y'=1:10), aes(x, y)) + geom_point()

#### Tests ---------------------------------------------------------------------
test_that('default theme', code={
  expect_equal_to_reference(theme_camprot(), 'theme_camprot.rds')
})

test_that('base_size = 10', code={
  expect_equal_to_reference(theme_camprot(base_size=10), 'theme_camprot_10.rds')
})

test_that('sans serif', code={
  expect_equal_to_reference(theme_camprot(base_family='sans'), 'theme_camprot_sans.rds')
})

test_that('not square', code={
  expect_equal_to_reference(theme_camprot(aspect_square=FALSE), 'theme_camprot_not_square.rds')
})

test_that('no border', code={
  expect_equal_to_reference(theme_camprot(border=FALSE), 'theme_camprot_no_border.rds')
})

#### Sanity checks -------------------------------------------------------------
