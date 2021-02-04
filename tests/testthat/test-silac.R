context('silac')

#### Setup ---------------------------------------------------------------------

#### Tests ---------------------------------------------------------------------
test_that("No heavy produces ratio = 0", {
  expect_equal(get_incorporation(light = 1, heavy = NA), 0)
})

test_that("No light produces ratio = 1", {
  expect_equal(get_incorporation(light = NA, heavy = 1), 1)
})

test_that("No heavy or light produces ratio = NA", {
  expect_equal(get_incorporation(light = NA, heavy = NA), NA)
})

test_that("Both heavy and light produces the correct ratio", {
  expect_equal(get_incorporation(light = 1, heavy = 2), 2/3)
})

#### Sanity checks -------------------------------------------------------------
