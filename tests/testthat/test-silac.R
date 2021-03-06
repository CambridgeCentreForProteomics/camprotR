context('silac')

#### Setup ---------------------------------------------------------------------
input <- data.frame('incorporation'=c(1, seq(.9,1,.01)))
output <- c(median(100*input$incorporation),
            mean(100*input$incorporation))

input_mix <- data.frame('incorporation'=seq(.45,.55,.01),
                        'corrected_incorporation'=seq(.95,1.05,.01))
output_mix <- c(median(100*input_mix$incorporation),
                median(100*input_mix$corrected_incorporation))

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

test_that("No mix", {
  expect_equal(unname(unlist(plot_incorporation(input)$'incorporation_estimates')),
               output)
})

test_that("Mix", {
  expect_equal(unname(unlist(
    plot_incorporation(input_mix, mix=.5)$'incorporation_estimates')),
    output_mix)
})

#### Sanity checks -------------------------------------------------------------
