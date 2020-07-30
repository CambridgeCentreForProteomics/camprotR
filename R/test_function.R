#' Test function
#'
#' Generates a scatter plot with ggplot2.
#'
#' @param data \code{data.frame}. Must have 2 numeric columns to plot.
#' @param x \code{numeric}. X values to plot.
#' @param y \code{numeric}. Y values to plot.
#'
#' @return ggplot using geom_point.
#' @example inst/examples/ex-test_function.R
#' @export
test_function <- function(data, x, y) {
  ggplot(data, aes(x = x, y = y)) +
    geom_point()
}
