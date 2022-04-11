#' Compute the indicator function
#'
#' Given a numeric vector, this function returns the logical vector indicating
#' whether each component is nonnegative or not.
#'
#' @param x A numeric vector.
#' @examples
#' compute_indicator(c(-1.6, 2.3, 0))
#' @export
compute_indicator <- function(x) {
  (x >= 0)
}


#' Compute the hinge function
#'
#' Given a numeric vector, this function evaluates the hinge function
#' \eqn{x -> max(x, 0)} for each component.
#'
#' @param x A numeric vector.
#' @examples
#' compute_hinge(c(-1.6, 2.3, 0))
#' @export
compute_hinge <- function(x) {
  pmax(x, 0)
}
