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
#' @param is_lattice A logical scalar for whether the design is lattice or not.
#'   Only used for "em", "hk", and "emhk".
#' @param number_of_bins An integer or an integer vector of the numbers of bins
#'   for the approximate methods. An integer if the numbers of bins are the same
#'   for all covariates. `NULL` if the approximate methods are not used.
#'   Currently available for "tc", "mars", and "tcmars".
#' @param extra_linear_covariates An integer vector or a string vector of extra
#'   linear covariates added to the model. Possibly used for "tc", "mars", and
#'   "tcmars".
#' @param coefficients A numeric vector of the coefficients of basis functions
#'   in the fitted model (in a scaled domain). A vector obtained by removing the
#'   zero components from the solution to the LASSO problem.
#' @param is_nonzero_component A logical vector indicating whether or not each
#'   component of the solution to the LASSO problem is nonzero.
#' @param is_included_basis A logical vector indicating whether or not each
#'   basis function is included in the LASSO problem.
#' @references Ki, D., Fang, B., and Guntuboyina, A. (2024+). MARS via LASSO.
#'   Accepted at \emph{Annals of Statistics}. Available at
#'   \url{https://arxiv.org/abs/2111.11694}.
#' @references Fang, B., Guntuboyina, A., and Sen, B. (2021). Multivariate
#'   extensions of isotonic regression and total variation denoising via entire
#'   monotonicity and Hardy—Krause variation. \emph{Annals of Statistics},
#'   \strong{49}(2), 769-792.
compute_fit <- function(X_eval, X_design, s, method, is_lattice, number_of_bins,
                        extra_linear_covariates, coefficients,
                        is_nonzero_component, is_included_basis = NULL) {
  M <- get_lasso_matrix(X_eval, X_design, s, method, is_lattice, number_of_bins,
                        extra_linear_covariates, is_included_basis)$lasso_matrix

  if (!is.null(is_nonzero_component)) {
    if (length(is_nonzero_component) != ncol(M)) {
      stop('`length(is_nonzero_component)` must be equal to the number of columns of the LASSO matrix.')
    }

    M <- M[, is_nonzero_component, drop = FALSE]
  }

  M %*% coefficients
}
