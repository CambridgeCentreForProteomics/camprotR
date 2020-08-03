#' Remove duplicated full stops
#'
#' @description A convenience function to remove any duplicated full stops
#' (aka periods) from the elements of a vector, e.g. will change `..` or `...`
#' or `....` etc. to just `.` Mainly used to fix column names.
#'
#' @param x `character`. Contains duplicate full stops to be removed.
#'
#' @return Returns `character` with duplicate full stops removed.
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

#' @noRd
message_parse <- function(x, column, message) {
  message(sprintf("%s features found from %s master proteins => %s",
                  nrow(x), length(unique(x[[column]])), message))
}
