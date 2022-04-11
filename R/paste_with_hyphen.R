#' Paste two strings with a hyphen
#'
#' Given two strings of positive lengths, this function returns the
#' concatenation of two strings separated by a hyphen. If one of the given two
#' strings has zero length, this function returns the concatenation without
#' separation.
#'
#' @param a A numeric vector. An observation vector.
#' @param b A numeric vector. An estimate vector of the same length.
paste_with_hyphen <- function(a, b) {
  if (a == "" || b == "") {
    paste0(a, b)
  } else {
    paste(a, b, sep = "-")
  }
}
