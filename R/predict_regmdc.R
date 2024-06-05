#' Give predictions based on a nonparametric regression model with mixed derivative constraints
#'
#' Given the model built up from \code{\link{regmdc}}, this function provides
#' predictions at new data points.
#'
#' @param regmdc_model An object of the class "regmdc". It is an object returned
#'   by \code{\link{regmdc}}.
#' @param X_pred A numeric matrix of new data points. Each row corresponds to
#'   an individual data point at which prediction will be made.
#'
#' @seealso \code{\link{regmdc}}, which produces nonparametric regression models
#'   with mixed derivative constraints fit to data.
#' @examples
#' fstar <- function(x) {x[1]**2 + x[2]**2}
#' X_design <- expand.grid(rep(list(seq(0, 13.0/14, length.out = 14L)), 2L))
#' theta <- apply(X_design, MARGIN = 1L, FUN = fstar)
#' sigma <- 1.0
#' y <- theta + sigma * rnorm(nrow(X_design))
#'
#' mars_model <- regmdc(X_design, y, s = 2L, method = "mars", V = 4.0, is_lattice = TRUE)
#'
#' X_pred <- matrix(c(1.0/3, 2.0/3, 2.0/3, 1.0/3), nrow = 2L)
#' predict_regmdc(mars_model, X_pred)
#' @export
predict_regmdc <- function(regmdc_model, X_pred) {
  # Error handling =============================================================
  if (!inherits(regmdc_model, "regmdc")) {
    stop('`regmdc_model` must be an object of the class "regmdc".')
  }
  # ============================================================================

  X_design <- regmdc_model$X_design
  s <- regmdc_model$s
  method <- regmdc_model$method
  is_lattice <- regmdc_model$is_lattice
  number_of_bins <- regmdc_model$number_of_bins
  extra_linear_covariates <- regmdc_model$extra_linear_covariates
  compressed_solution <- regmdc_model$compressed_solution
  is_nonzero_component <- regmdc_model$is_nonzero_component
  is_included_basis <- regmdc_model$is_included_basis

  predicted_values <- compute_fit(X_pred, X_design, s, method, is_lattice,
                                  number_of_bins, extra_linear_covariates,
                                  compressed_solution, is_nonzero_component,
                                  is_included_basis)

  predicted_values
}
