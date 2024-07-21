#' Give predictions based on a nonparametric regression model with mixed derivative constraints
#'
#' Given the model built up from \code{\link{regmdc}}, this function provides
#' predictions at new data points.
#'
#' @param regmdc_model An object of the class "regmdc". It is an object returned
#'   by \code{\link{regmdc}}.
#' @param X_pred A numeric matrix of new data points. Each row corresponds to
#'   an individual data point at which prediction will be made.
#' @references Ki, D., Fang, B., and Guntuboyina, A. (2024+). MARS via LASSO.
#'   Accepted at \emph{Annals of Statistics}. Available at
#'   \url{https://arxiv.org/abs/2111.11694}.
#' @references Fang, B., Guntuboyina, A., and Sen, B. (2021). Multivariate
#'   extensions of isotonic regression and total variation denoising via entire
#'   monotonicity and Hardy—Krause variation. \emph{Annals of Statistics},
#'   \strong{49}(2), 769-792.
#' @seealso \code{\link{regmdc}}, which produces nonparametric regression models
#'   with mixed derivative constraints fit to data.
#' @examples
#' fstar <- function(x) {(
#'   (x[1] - 0.25 >= 0) + (x[2] - 0.25 >= 0)
#'   + (x[1] - 0.25 >= 0) * (x[2] - 0.25 >= 0)
#' )}  # the true underlying function
#' X_design <- expand.grid(rep(list(seq(0, 1, length.out = 10L)), 3L))
#' theta <- apply(X_design, MARGIN = 1L, FUN = fstar)
#' sigma <- 0.1
#' y <- theta + sigma * rnorm(nrow(X_design))
#'
#' hk_model <- regmdc(X_design, y, s = 2L, method = "hk", V = 3.0)
#'
#' X_pred <- c(1.0/3, 2.0/3, 1.0/3)
#' predict_regmdc(hk_model, X_pred)
#' X_pred <- matrix(c(1.0/3, 2.0/3, 1.0/3,
#'                    2.0/3, 1.0/3, 2.0/3),
#'                  ncol = 3L, byrow = TRUE)
#' predict_regmdc(hk_model, X_pred)
#'
#' fstar <- function(x) {(
#'   - max(x[1] - 0.25, 0) - max(x[2] - 0.25, 0)
#'   - max(x[1] - 0.25, 0) * max(x[2] - 0.25, 0)
#' )}  # the true underlying function
#' X_design <- expand.grid(rep(list(seq(0, 1, length.out = 10L)), 3L))
#' theta <- apply(X_design, MARGIN = 1L, FUN = fstar)
#' sigma <- 0.1
#' y <- theta + sigma * rnorm(nrow(X_design))
#'
#' mars_model <- regmdc(X_design, y, s = 2L, method = "mars", V = 3.0)
#'
#' X_pred <- c(1.0/3, 2.0/3, 1.0/3)
#' predict_regmdc(mars_model, X_pred)
#' X_pred <- matrix(c(1.0/3, 2.0/3, 1.0/3,
#'                    2.0/3, 1.0/3, 2.0/3),
#'                  ncol = 3L, byrow = TRUE)
#' predict_regmdc(mars_model, X_pred)
#' @export
predict_regmdc <- function(regmdc_model, X_pred) {
  X_design <- regmdc_model$X_design
  max_vals <- regmdc_model$max_vals
  min_vals <- regmdc_model$min_vals
  s <- regmdc_model$s
  method <- regmdc_model$method
  is_scaled <- regmdc_model$is_scaled
  is_lattice <- regmdc_model$is_lattice
  number_of_bins <- regmdc_model$number_of_bins
  increasing_covariates <- regmdc_model$increasing_covariates
  decreasing_covariates <- regmdc_model$decreasing_covariates
  concave_covariates <- regmdc_model$concave_covariates
  convex_covariates <- regmdc_model$convex_covariates
  variation_constrained_covariates <- regmdc_model$variation_constrained_covariates
  extra_linear_covariates <- regmdc_model$extra_linear_covariates
  is_included_basis <- regmdc_model$is_included_basis
  is_nonzero_component <- regmdc_model$is_nonzero_component
  coefficients <- regmdc_model$coefficients

  # Error handling =============================================================
  if (!inherits(regmdc_model, "regmdc")) {
    stop('`regmdc_model` must be an object of the class "regmdc".')
  }

  if (is.vector(X_pred)) {
    if (length(X_pred) != ncol(X_design)) {
      stop('`length(X_pred)` must be equal to `ncol(X_design)`.')
    }

    X_pred <- matrix(X_pred, nrow = 1)
  }
  if (anyNA(X_pred)) {
    stop('`X_pred` must not include NA or NaN.')
  }
  if (!all(apply(X_pred, MARGIN = c(1L, 2L), FUN = is.numeric))) {
    stop('`X_pred` must be numeric.')
  }
  if (ncol(X_pred) != ncol(X_design)) {
    stop('`ncol(X_pred)` must be equal to `ncol(X_design)`.')
  }

  # ============================================================================
  M <- get_lasso_matrix(
    X_pred, X_design, max_vals, min_vals, s, method, is_scaled, is_lattice,
    number_of_bins, increasing_covariates, decreasing_covariates,
    concave_covariates, convex_covariates, variation_constrained_covariates,
    extra_linear_covariates, is_included_basis
  )$lasso_matrix

  if (!is.null(is_nonzero_component)) {
    if (length(is_nonzero_component) != ncol(M)) {
      stop(strwrap('`length(is_nonzero_component)` must be equal to the number
                   of columns of the LASSO matrix.'))
    }

    M <- M[, is_nonzero_component, drop = FALSE]
  }

  M %*% coefficients
}
