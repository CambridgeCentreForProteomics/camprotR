#' Get ratio of two numbers
#'
#' Calculates the ratio of two `numeric` columns
#' in a data.frame, and indicates if there are missing values. Calculation is:
#' `trt - ctr`. Typically used to calculate the ratio of two `log2`
#' transformed SILAC/label-free intensity columns.
#'
#' @param data `data.frame`.
#' @param treatment,control `variable`. Name of numeric column in the data.frame.
#' @param bind `logical`. Should the resulting ratios be added to data
#' as a column? Default is `TRUE`.
#'
#' @return Returns `data` with 2 new columns: ratio (a numeric column) and
#' missing (a character column indicating if there was missing values in
#' `treament` and/or `control`). Returns a `named list` if `bind = FALSE`.
#' @examples
#' my_data <- data.frame(
#'   treatment = c(30, NA, 12, NA),
#'   control = c(20, 12, NA, NA)
#' )
#'
#' my_data2 <- get_ratio(my_data, treatment, control)
#' @export
get_ratio <- function(data, treatment, control, bind = TRUE) {
  treatment <- deparse(substitute(treatment))
  control <- deparse(substitute(control))

  ratios <- apply(
    data[, c(treatment, control)],
    MARGIN = 1,
    FUN = function(x) {
      ratio <- x[1] - x[2]
      missing <- sum(is.na(x))

      if (missing == 1) {
        if (is.na(x[1])) missing <- sprintf("%s Missing", names(x[1]))
        else if (is.na(x[2])) missing <- sprintf("%s Missing", names(x[2]))
      } else if (missing == 0) {
        missing <- "Neither missing"
      } else if (missing == 2) {
        missing <- "Both missing"
      }

      c(ratio, missing)
    }
  )
  result <- list(
    ratio = as.numeric(ratios[1, ]),
    missing = as.factor(ratios[2, ])
  )
  if (bind) result <- cbind(data, result)
  result
}
