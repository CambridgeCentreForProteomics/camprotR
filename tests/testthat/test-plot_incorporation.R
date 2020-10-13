context('Plot incorporation')

# set up #
input <- data.frame('incorporation'=c(1, seq(.9,1,.01)))
output <- c(median(100*input$incorporation),
            mean(100*input$incorporation))

input_mix <- data.frame('incorporation'=seq(.45,.55,.01),
                        'corrected_incorporation'=seq(.95,1.05,.01))
output_mix <- c(median(100*input_mix$incorporation),
                median(100*input_mix$corrected_incorporation))
##########

# tests #
test_that("No mix", {
  expect_equal(unname(unlist(plot_incorporation(input)$'incorporation_estimates')),
               output)
})

# tests #
test_that("Mix", {
  expect_equal(unname(unlist(
    plot_incorporation(input_mix, mix=.5)$'incorporation_estimates')),
    output_mix)
})
