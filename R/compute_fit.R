#' Compute the fitted values of an estimation method
#'
#' Given the solution to the corresponding LASSO problem obtained from
#' \code{\link{get_lasso_problem_soln}}, this function computes
#' the fitted values of an estimation method. Available estimation methods are
#' entirely monotonic regression ("em"), Hardy—Krause variation denoising ("hk"),
#' their generalization ("emhk"), totally convex regression ("tc"), MARS via
#' LASSO ("mars"), and their generalization ("tcmars").
#'
#' @param X_eval A numeric evaluation matrix. Each row corresponds to individual
#'   evaluation point at which the fitted value of the estimation method is
#'   computed.
#' @param X_design A numeric design matrix. Each row corresponds to individual
#'   data.
#' @param s A numeric scalar indicating the maximum order of interaction between
#'   covariates allowed in the estimation method.
#' @param method A string indicating the estimation method. One of "em", "hk",
#'   "emhk", "tc", "mars", and "tcmars".
#' @param compressed_soln A numeric vector obtained by removing the zero
#'   components from the solution to the LASSO problem.
#' @param is_nonzero A logical vector indicating whether each component of the
#'   solution to the LASSO problem is nonzero or not.
#' @details
#' You can use this function without removing the zero components of the
#' original solution. You can simply put the solution to the LASSO problem to
#' `compressed_soln` and ignore `is_nonzero` (or set `is_nonzero` = NULL).
#' @seealso \code{\link{get_lasso_problem_soln}}, which is used for finding the
#'   solution to the LASSO problem.
#' @references Ki, D., Fang, B., and Guntuboyina, A. (2021). MARS via LASSO.
#'   \url{https://arxiv.org/abs/2111.11694}.
#' @references Fang, B., Guntuboyina, A., and Sen, B. (2021). Multivariate
#' extensions of isotonic regression and total variation denoising via entire
#' monotonicity and Hardy—Krause variation. \emph{The Annals of Statistics},
#' \strong{49}(2), 769-792.
#' @examples
#' fstar <- function(x) {x[1] + x[2]}
#' X_design <- expand.grid(rep(list(seq(0, 9.0/10, length.out = 10L)), 2))
#' theta <- apply(X_design, MARGIN = 1L, FUN = fstar)
#' sigma <- 1
#' y <- theta + sigma * rnorm(nrow(X_design))
#' raw_soln <- get_lasso_problem_soln(X_design, y, s = 1L, method = "em")$solution
#' X_eval <- expand.grid(rep(list(seq(0, 1, length.out = 51L)), 2L))
#'
#' compute_fit(X_eval, X_design, s = 1L, method = "em", compressed_soln = raw_soln)
#'
#' threshold <- 1e-4
#' is_nonzero <- (abs(raw_soln) > threshold)
#' processed_soln <- as.numeric(is_nonzero) * raw_soln
#' compressed_soln <- processed_soln[processed_soln != 0]
#' compute_fit(X_eval, X_design, s = 1L, method = "em", compressed_soln, is_nonzero)
#' @export
compute_fit <- function(X_eval, X_design, s, method, compressed_soln,
                        is_nonzero = NULL) {
  # Error handling =============================================================
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

  if (!is.integer(s)) {
    stop('`s` must be an integer.')
  }
  if (s < 1L || s > ncol(X_design)) {
    stop('`s` must be at least 1 and at most `ncol(X_design)`.')
  }

  if (!(method %in% c('em', 'hk', 'emhk', 'tc', 'mars', 'tcmars'))) {
    stop('`method` must be one of "em", "hk", "emhk", "tc", "mars", and "tcmars".')
  }

  if (!is.null(is_nonzero)) {
    if (!is.logical(is_nonzero)) {
      stop('`is_nonzero` must be logical.')
    }
    if (length(compressed_soln) != sum(as.numeric(is_nonzero))) {
      stop('`length(compressed_soln)` must be equal to the number of TRUEs in `is_nonzero`.')
    }
  }
  # ============================================================================

  M <- get_lasso_matrix(X_eval, X_design, s, method)$lasso_matrix

  if (!is.null(is_nonzero)) {
    if (length(is_nonzero) != ncol(M)) {
      stop('`length(is_nonzero)` must be equal to the number of columns of the LASSO matrix.')
    }

    M <- M[, is_nonzero]
  }

  M %*% compressed_soln
}
