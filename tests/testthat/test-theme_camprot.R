test_that("theme_camprot() works", {
  p <- ggplot(data.frame(x = c(1:10), y = c(1:10)), aes(x = x, y = y)) +
    geom_point()

  # test that defaults work
  vdiffr::expect_doppelganger(
    "theme-camprot-default",
    p + theme_camprot()
  )

  # test that you can change base size
  vdiffr::expect_doppelganger(
    "theme-camprot-base-size-10",
    p + theme_camprot(base_size = 10)
  )

  # test that you can change the font family
  vdiffr::expect_doppelganger(
    "theme-camprot-sans-serif",
    p + theme_camprot(base_family = "sans")
  )

  # test that you can turn off square aspect ratio
  vdiffr::expect_doppelganger(
    "theme-camprot-not-square",
    p + theme_camprot(aspect_square = FALSE)
  )

  # test that you can turn off the border
  vdiffr::expect_doppelganger(
    "theme-camprot-no-border",
    p + theme_camprot(border = FALSE)
  )
})
