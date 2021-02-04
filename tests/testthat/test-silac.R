context('SILAC incorporation')

test_that("No heavy", {
  expect_equal(get_incorporation(light = 1, heavy = NA), 0)
})

test_that("No light", {
  expect_equal(get_incorporation(light = NA, heavy = 1), 1)
})

test_that("No heavy or light", {
  expect_equal(get_incorporation(light = NA, heavy = NA), NA)
})

test_that("Both heavy and light", {
  expect_equal(get_incorporation(light = 1, heavy = 2), 2/3)
})
