context('get_cat_palette')

#### Setup ---------------------------------------------------------------------

#### Tests ---------------------------------------------------------------------
test_that("zero colours", {
  expect_error(get_cat_palette(0))
})

test_that("13 colours", {
  expect_error(get_cat_palette(13))
})

test_that("1 colour", {
  expect_equal(get_cat_palette(1), "#2271B2")
})

test_that("7 colours", {
  expect_equal(get_cat_palette(7),
               c("#2271B2", "#d55e00", "#359B73", "#e69f00",
                 "#3DB7E9", "#f0e442", "#F748A5"))
})


test_that("8 colours", {
  expect_equal(get_cat_palette(8),
               c("#9F0162", "#009F81", "#FF5AAF", "#00FCCF",
                 "#8400CD", "#008DF9", "#00C2F9", "#FFB2FD"))
})


test_that("12 colours", {
  expect_equal(get_cat_palette(12),
               c("#9F0162", "#009F81", "#FF5AAF", "#00FCCF",
                 "#8400CD", "#008DF9", "#00C2F9", "#FFB2FD",
                 "#A40122", "#E20134", "#FF6E3A", "#FFC33B"))
})

#### Sanity checks -------------------------------------------------------------
