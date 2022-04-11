#' Compute the mean squared error of an estimator
#'
#' @param y A numeric observation vector.
#' @param y_hat A numeric estimate vector of the same length.
#' @examples
#' compute_mse(c(1.0, 2.3, 0.8), c(1.3, 1.9, 0.7))
#' @export
compute_mse <- function(y, y_hat) {
  # Error handling =============================================================
  if (anyNA(y)) {
    stop('`y` must not include NA or NaN')
  }
  if (anyNA(y_hat)) {
    stop('`y_hat` must not include NA or NaN')
  }
  if (!is.numeric(y)) {
    stop('`y` must be a numeric vector')
  }
  if (!is.numeric(y_hat)) {
    stop('`y_hat` must be a numeric vector')
  }
  if (length(y) != length(y_hat)) {
    stop('`y` and `y_hat` must have the same length')
  }
  # ============================================================================

  mean((y - y_hat)**2)
}
