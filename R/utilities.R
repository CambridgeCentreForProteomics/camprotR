#' Remove duplicated full stops
#'
#' @description A convenience function to remove any duplicated full stops
#' (aka periods) from the elements of a vector, e.g. will change `..` or `...`
#' or `....` etc. to just `.` Mainly used to fix column names.
#'
#' @param x `character` or `string`. Contains duplicate full stops to be removed.
#'
#' @return Returns `character` or `string` with duplicate full stops removed.
#' @examples
#'
#' df <- data.frame(
#'   column...name = c(1, 2, 3)
#' )
#'
#' colnames(df) <- remove_dots(colnames(df))
#'
#' @export
remove_dots <- function(x) {
  gsub("(?<=\\.)\\.+", "", x, perl = TRUE)
}

#' Remove leading X
#'
#' @description A convenience function to remove a leading capital X.
#' Is case sensitive.
#'
#' @param x `character` or `string`.
#' @return Returns `character` or `string` with leading X removed.
#' @examples
#'
#' df <- data.frame('X1'=c(1,2))
#'
#' removeX(colnames(df))
#'
#' @export
remove_x <- function(x) {
  gsub("^X", "", x)
}

#' @noRd
message_parse <- function(x, column, message) {
  message(sprintf("%s features found from %s master proteins => %s",
                  nrow(x), length(unique(x[[column]])), message))
}
