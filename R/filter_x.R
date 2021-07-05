#' @noRd
subset_na <- function(data, op = c("==", "<=", ">=", "!="), regex, value) {
  # refer to operator by name
  op <- as.name(op)

  # get number of NA in column group
  n_na <- apply(
    is.na(data[, c(grep(regex, colnames(data)))]),
    1,
    sum
  )

  # subset data based on number of NA in column group
  subset(data, sapply(n_na, op, value))
}

#' Filter rows based on `NA` values in column groups
#'
#' @description This function is used to filter rows based on number of `NA`
#' values in different groups of columns. The column groups are defined by a
#' regular expression matching their names. This provides more granular control
#' over what is filtered out compared to filtering based on ratios alone.
#'
#' The function can be explained as follows:
#' ```
#' filter_na(
#'   data,
#'   op = "<=",
#'   setop = intersect,
#'   regex = c("sample_A", "sample_B"),
#'   value = c(1, 1)
#' )
#' ```
#'
#' The above can be translated as: keep rows that have less
#' than or equal to 1 `NA` value in the "sample_A" columns AND the "sample_B"
#' columns.
#'
#' More illustrative examples are below.
#'
#' @param data `data.frame` or `matrix` to be filtered.
#' @param op `string`. Defines the operation used for the filtering. Should be
#' one of `"=="`, `"<="`, `">="`, or `"!="`.
#' @param setop The function used for combining the results for the different
#' column groups. Should be one of `union` (equivalent to "OR"), `intersect`
#' (equivalent to "AND"), or `setdiff` (equivalent to "SYMMETRIC DIFFERENCE").
#' @param regex `character vector` of length n. The regular expression(s) used
#' to define the column groups. The length of the vector indicates the number
#' of column groups to use.
#' @param value `numeric vector` of length n. The number of `NA` values to check
#' for in the rows of each column group. Must be the same length as the
#' `regex` vector.
#'
#' @return Returns a filtered object of the same class as `data`.
#' @examples
#' library(magrittr)
#'
#' # set a seed for reproducibility
#' set.seed(123)
#'
#' # make a data.frame with fake protein/peptide abundance data
#' df <- cbind(replicate(3, runif(7, 13, 24)),
#'             replicate(3, runif(7, 24, 30))) %>%
#'  as.data.frame() %>%
#'  `colnames<-`(c(paste0("ctr", 1:3), paste0("trt", 1:3)))
#'
#' # add in some missing values
#' mapply(function(i, j) {
#'   df[i, j] <<- NA_real_
#' }, c(2:7), list(2, c(2:3), c(3:4), c(2:5), c(2:6), c(1:6)))
#'
#' # filter for rows with <= 1 NA value in ctr OR trt samples
#' filter_na(df, "<=", union, c("ctr", "trt"), c(1, 1))
#'
#' # filter for rows with <= 1 NA value in ctr AND trt samples
#' filter_na(df, "<=", intersect, c("ctr", "trt"), c(1, 1))
#'
#' # filter for rows with exactly 1 NA value in ctr AND trt samples
#' filter_na(df, "==", intersect, c("ctr", "trt"), c(1, 1))
#'
#' # filter for rows with at least 1 non-NA value
#' filter_na(df, "<=", union, "ctr|trt", 5)
#'
#' # filter for rows with at least 3 NA values
#' filter_na(df, ">=", union, "ctr|trt", 3)
#'
#' @export
filter_na <- function(data, op, setop, regex, value) {
  # repeat subset_na for pairs of values in two vectors
  # then combine the result depending on the desired set operation
  result <- mapply(function(r, v) subset_na(data, op, r, v),
                   regex, value, SIMPLIFY = FALSE)
  Reduce(setop, result)
}

#' @noRd
subset_zero <- function(data, op = c("==", "<=", ">=", "!="), regex, value) {
  # refer to operator by name
  op <- as.name(op)

  # get sum of zero values in column group
  n_zero <- apply(
    data[, c(grep(regex, colnames(data)))],
    1,
    function(x) sum(x == 0, na.rm = TRUE)
  )

  # subset data based on number of zeros in column group
  subset(data, sapply(n_zero, op, value))
}

#' Filter rows based on zero values in column groups
#'
#' @description This function is used to filter rows based on number of zero
#' values in different groups of columns. The column groups are defined by a
#' regex matching their names. This provides more granular control over what
#' is filtered out compared to filtering based on ratios alone.
#'
#' The function can be explained as follows:
#' ```
#' filter_zero(
#'   data,
#'   op = "<=",
#'   setop = intersect,
#'   regex = c("sample_A", "sample_B"),
#'   value = c(1, 1)
#' )
#' ```
#'
#' The above can be translated as: keep rows that have less
#' than or equal to 1 zero value in the "sample_A" columns AND the "sample_B"
#' columns.
#'
#' More illustrative examples are below.
#'
#' @param data `data.frame` or `matrix` to be filtered.
#' @param op `string`. Defines the operation used for the filtering. Should be
#' one of `"=="`, `"<="`, `">="`, or `"!="`.
#' @param setop The function used for combining the results for the different
#' column groups. Should be one of `union` (equivalent to "OR"), `intersect`
#' (equivalent to "AND"), or `setdiff` (equivalent to "SYMMETRIC DIFFERENCE").
#' @param regex `character vector` of length n. The regular expression(s) used
#' to define the column groups. The length of the vector indicates the number
#' of column groups to use.
#' @param value `numeric vector` of length n. The number of zero values to check
#' for in the rows of each column group. Must be the same length as the
#' `regex` vector.
#'
#' @return Returns a filtered object of the same class as `data`.
#' @examples
#' #
#' @export
filter_zero <- function(data, op, setop, regex, value) {
  # repeat subset_zero for pairs of values in two vectors
  # then combine the result depending on the desired set operation
  result <- mapply(function(r, v) subset_zero(data, op, r, v),
                   regex, value, SIMPLIFY = FALSE)
  Reduce(setop, result)
}

#' @noRd
subset_val <- function(data, op = c("==", "<=", ">=", "!="), regex, value) {
  # refer to operator by name
  op <- as.name(op)

  # get sum of values in column group
  col_sum <- apply(
    data[, c(grep(regex, colnames(data)))],
    1,
    function(x) sum(x, na.rm = TRUE)
  )

  # subset data based on sum of values in column group
  subset(data, sapply(col_sum, op, value))
}

#' Filter rows based on the sum of values in column groups
#'
#' @description This function is used to filter rows based on the sum of
#' values in different groups of columns. The column groups are defined by a
#' regex matching their names. This provides more granular control over what
#' is filtered out compared to filtering based on ratios alone.
#'
#' The function can be explained as follows:
#' ```
#' filter_val(
#'   data,
#'   op = "<=",
#'   setop = intersect,
#'   regex = c("sample_A", "sample_B"),
#'   value = c(100, 50)
#' )
#' ```
#'
#' The above can be translated as: keep rows that add up to less than or equal
#' to 100 in the "sample_A" columns AND less than or equal to 50 in the
#' "sample_B" columns.
#'
#' More illustrative examples are below.
#'
#' @param data `data.frame` or `matrix` to be filtered.
#' @param op `string`. Defines the operation used for the filtering. Should be
#' one of `"=="`, `"<="`, `">="`, or `"!="`.
#' @param setop The function used for combining the results for the different
#' column groups. Should be one of `union` (equivalent to "OR"), `intersect`
#' (equivalent to "AND"), or `setdiff` (equivalent to "SYMMETRIC DIFFERENCE").
#' @param regex `character vector` of length n. The regular expression(s) used
#' to define the column groups. The length of the vector indicates the number
#' of column groups to use.
#' @param value `numeric vector` of length n. The number of summed values to
#' check for in the rows of each column group. Must be the same length as the
#' `regex` vector.
#'
#' @return Returns a filtered object of the same class as `data`.
#' @examples
#' #
#' @export
filter_val <- function(data, op, setop, regex, value) {
  # repeat subset_val for pairs of values in two vectors
  # then combine the result depending on the desired set operation
  result <- mapply(function(r, v) subset_val(data, op, r, v),
                   regex, value, SIMPLIFY = FALSE)
  Reduce(setop, result)
}
