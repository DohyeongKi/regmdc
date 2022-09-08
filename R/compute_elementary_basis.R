#' Compute the indicator function
#'
#' Given a numeric vector, this function returns the logical vector indicating
#' whether or not each component is nonnegative.
#'
#' @param x A numeric vector.
compute_indicator <- function(x) {
  (x >= 0)
}


#' Compute the hinge function
#'
#' Given a numeric vector, this function evaluates the hinge function
#' \eqn{x -> max(x, 0)} for each component.
#'
#' @param x A numeric vector.
compute_hinge <- function(x) {
  pmax(x, 0)
}
