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

test_that("plot_incorporation() works with no mixing", {
  fake_HL_ratios <- data.frame(
    incorporation = c(1, seq(0.9, 1, 0.01))
  )

  # target is a list of length == 2
  # list contains a plot and another list
  target <- plot_incorporation(
    data = fake_HL_ratios,
    incorporation_col = "incorporation",
    level = "peptide",
    mix = 0
  )

  # test that incorporation plot output works
  vdiffr::expect_doppelganger(
    "plot-incorporation-no-mix",
    plot(target$p)
  )

  # test that list output is okay
  expect_equal(
    target$incorporation_estimates,
    list(
      "Median incorporation" = median(fake_HL_ratios$incorporation * 100),
      "Mean incorporation" = mean(fake_HL_ratios$incorporation * 100)
    )
  )
})

test_that("plot_incorporation() works with mixing", {
  fake_HL_ratios <- data.frame(
    incorporation = seq(0.45, 0.55, 0.01),
    corrected_incorporation = seq(0.95, 1.05, 0.01)
  )

  # target is a list of length == 2
  # list contains a plot and another list
  target <- plot_incorporation(
    data = fake_HL_ratios,
    incorporation_col = "incorporation",
    corrected_col = "corrected_incorporation",
    level = "peptide",
    mix = 0.5
  )

  # test that incorporation plot output works
  vdiffr::expect_doppelganger(
    "plot-incorporation-mix",
    plot(target$p)
  )

  # test that list output is okay
  expect_equal(
    target$incorporation_estimates,
    list(
      "Median incorporation" = median(fake_HL_ratios$incorporation * 100),
      "Median incorporation (corrected)" = median(fake_HL_ratios$corrected_incorporation * 100)
    )
  )
})
