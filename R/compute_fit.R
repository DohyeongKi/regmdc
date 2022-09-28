#' Compute the fitted values of an estimation method at given points.
#'
#' Given the solution to the corresponding LASSO problem, this function computes
#' the fitted values of an estimation method.
#'
#' @param X_eval A numeric evaluation matrix. Each row corresponds to an
#'   individual evaluation point at which the fitted value of the estimation
#'   method is computed.
#' @param X_design A numeric design matrix. Each row corresponds to an
#'   individual data.
#' @param s A numeric scalar indicating the maximum order of interaction between
#'   covariates allowed in the estimation method.
#' @param method A string indicating the estimation method. One of "em", "hk",
#'   "emhk", "tc", "mars", and "tcmars".
#' @param number_of_bins An integer or an integer vector of the numbers of bins
#'   for the approximate methods. An integer if the numbers of bins are the same
#'   for all covariates. `NULL` if the approximate methods are not used.
#' @param compressed_solution A numeric vector obtained by removing the zero
#'   components from the solution to the LASSO problem.
#' @param is_nonzero_component A logical vector indicating whether or not each
#'   component of the solution to the LASSO problem is nonzero.
#' @param is_included_basis A logical vector indicating whether or not each
#'   basis function is included in the LASSO problem.
#' @references Ki, D., Fang, B., and Guntuboyina, A. (2021). MARS via LASSO.
#'   \url{https://arxiv.org/abs/2111.11694}.
#' @references Fang, B., Guntuboyina, A., and Sen, B. (2021). Multivariate
#'   extensions of isotonic regression and total variation denoising via entire
#'   monotonicity and Hardy—Krause variation. \emph{The Annals of Statistics},
#'   \strong{49}(2), 769-792.
compute_fit <- function(X_eval, X_design, s, method, number_of_bins,
                        compressed_solution, is_nonzero_component,
                        is_included_basis = NULL) {
  # Error handling =============================================================
  if (length(dim(X_eval)) != 2L) {
    stop('`X_eval` must be a matrix or a data frame.')
  }
  if (length(dim(X_design)) != 2L) {
    stop('`X_design` must be a matrix or a data frame.')
  }
  if (anyNA(X_eval)) {
    stop('`X_eval` must not include NA or NaN.')
  }
  if (anyNA(X_design)) {
    stop('`X_design` must not include NA or NaN.')
  }
  if (!all(apply(X_eval, MARGIN = c(1L, 2L), FUN = is.numeric))) {
    stop('`X_eval` must be numeric.')
  }
  if (!all(apply(X_design, MARGIN = c(1L, 2L), FUN = is.numeric))) {
    stop('`X_design` must be numeric.')
  }
  if (min(X_design) < 0 || max(X_design) > 1) {
    stop('Every component of `X_design` must be between 0 and 1.')
  }
  if (ncol(X_eval) != ncol(X_design)) {
    stop('`ncol(X_eval)` must be equal to `ncol(X_design)`.')
  }

  if (!is.integer(s) || s > ncol(X_design)) {
    stop('`s` must be an integer.')
  }
  if (s < 1L || s > ncol(X_design)) {
    stop('`s` must be at least 1 and at most `ncol(X_design)`.')
  }

  if (!(method %in% c('em', 'hk', 'emhk', 'tc', 'mars', 'tcmars'))) {
    stop('`method` must be one of "em", "hk", "emhk", "tc", "mars", and "tcmars".')
  }

  if (!is.null(number_of_bins)) {
    if (!is.integer(number_of_bins)) {
      stop('`number_of_bins` must be an integer or an integer vector.')
    }

    if (length(number_of_bins) == 1L) {
      number_of_bins <- rep(number_of_bins, ncol(X_design))
    } else if (length(number_of_bins) != ncol(X_design)) {
      stop('`length(number_of_bins)` must be equal to 1 or `ncol(X_design)`.')
    }

    if (method %in% c('em', 'hk', 'emhk')) {
      stop('The approximate method is only available for `tc`, `mars`, and `tcmars` at this point.')
    }
  }

  if (!is.null(is_nonzero_component)) {
    if (!is.logical(is_nonzero_component)) {
      stop('`is_nonzero_component` must be logical.')
    }
    if (length(compressed_solution) != sum(as.numeric(is_nonzero_component))) {
      stop('`length(compressed_solution)` must be equal to the number of TRUEs in `is_nonzero_component`.')
    }
  }

  if (!is.null(is_included_basis) && !is.logical(is_included_basis)) {
    stop('`is_included_basis` must be logical.')
  }
  # ============================================================================

  M <- get_lasso_matrix(X_eval, X_design, s, method, number_of_bins,
                        is_included_basis)$lasso_matrix

  if (!is.null(is_nonzero_component)) {
    if (length(is_nonzero_component) != ncol(M)) {
      stop('`length(is_nonzero_component)` must be equal to the number of columns of the LASSO matrix.')
    }

    M <- M[, is_nonzero_component]
  }

  M %*% compressed_solution
}
