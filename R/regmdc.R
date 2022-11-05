#' Fit to data a nonparametric regression model with mixed derivative constraints
#'
#' Given an estimation method, this function builds up a model by solving the
#' corresponding constrained LASSO problem. Available estimation methods are
#' entirely monotonic regression ("em"), Hardy—Krause variation denoising ("hk"),
#' their generalization ("emhk"), totally convex regression ("tc"), MARS via
#' LASSO ("mars"), and their generalization ("tcmars"). For details on the
#' corresponding LASSO problems, see, for example, Section 3 of Fang et al.
#' (2021) (for entirely monotonic regression and Hardy—Krause variation
#' denoising) and Section 3 and 6 of Ki et al. (2021) (for MARS via LASSO).
#' There are some other ongoing research and working papers, and they would be
#' available in the future.
#'
#' @param X_design A numeric design matrix. Each row corresponds to an
#'   individual data.
#' @param y A numeric observation vector.
#' @param s A numeric scalar indicating the maximum order of interaction between
#'   covariates allowed in the estimation method.
#' @param method A string indicating the estimation method. One of "em", "hk",
#'   "emhk", "tc", "mars", and "tcmars".
#' @param V A numeric scalar. An upper bound on complexity measure (variation).
#'   Required for "hk" and "mars", and possibly used for "emhk" and "tcmars".
#' @param threshold A numeric scalar to determine whether each component of the
#'   solution to the LASSO problem is zero or not.
#' @param number_of_bins An integer or an integer vector of the numbers of bins
#'   for the approximate methods. An integer if the numbers of bins are the same
#'   for all covariates. `NULL` if the approximate methods are not used.
#' @param constrained_interactions A string vector indicating constrained
#'   interactions between covariates. Possibly used for "emhk" and "tcmars".
#' @param positive_interactions A string vector indicating positive interactions
#'   between covariates. Possibly used for "emhk".
#' @param negative_interactions A string vector indicating negative interactions
#'   between covariates. Possibly used for "emhk".
#' @param increasing_interactions A string vector indicating monotonically
#'   increasing interactions between covariates. Possibly used for "tcmars".
#' @param decreasing_interactions A string vector indicating monotonically
#'   decreasing interactions between covariates. Possibly used for "tcmars".
#' @details
#' Contrary to entirely monotonic regression (resp. totally convex regression)
#' where every interaction between covariates is restricted to be positive (resp.
#' monotonically increasing) and Hardy—Krause variation denoising (resp. MARS via
#' LASSO) where every interaction is constrained, in their generalization "emhk"
#' (resp. "tcmars"), you can freely choose which interaction to constrain and
#' which interaction to restrict to be positive (resp. monotonically increasing)
#' or negative (resp. monotonically decreasing). When you use any of the
#' interaction arguments, you need to input a vector of interactions each of
#' which is represented in a form of concatenation of covariate numbers, sorted
#' in an increasing order and separated by a hyphen. For example, if you want to
#' constrain interaction between the second, the third, and the fifth covariates,
#' you need to include "2-3-5" to your input vector of `constrained_interactions`.
#'
#' If `number_of_bins` is not `NULL`, then the approximate methods are used.
#' Refer to \code{\link{get_unique_column_entries}} to see what differences are
#' made in the approximate methods. The approximate methods for totally convex
#' regression, MARS via LASSO, and their generalization have been implemented.
#' Approximate methods for entirely monotonic regression, Hardy—Krause variation
#' denoising, and their generalization will be available in the future.
#'
#' @references Ki, D., Fang, B., and Guntuboyina, A. (2021). MARS via LASSO.
#'   \url{https://arxiv.org/abs/2111.11694}.
#' @references Fang, B., Guntuboyina, A., and Sen, B. (2021). Multivariate
#'   extensions of isotonic regression and total variation denoising via entire
#'   monotonicity and Hardy—Krause variation. \emph{The Annals of Statistics},
#'   \strong{49}(2), 769-792.
#' @examples
#' fstar <- function(x) {x[1]**2 + x[2]**2}
#' X_design <- expand.grid(rep(list(seq(0, 13.0/14, length.out = 14L)), 2L))
#' theta <- apply(X_design, MARGIN = 1L, FUN = fstar)
#' sigma <- 1.0
#' y <- theta + sigma * rnorm(nrow(X_design))
#'
#' regmdc(X_design, y, s = 1L, method = "em")
#' regmdc(X_design, y, s = 2L, method = "hk", V = 2.0)
#' regmdc(X_design, y, s = 2L, method = "emhk", V = 1.0,
#'        constrained_interactions = c('1-2'),
#'        positive_interactions = c('1'),
#'        negative_interactions = c('2'))
#' regmdc(X_design, y, s = 2L, method = "tc")
#' regmdc(X_design, y, s = 2L, method = "mars", V = 4.0)
#' regmdc(X_design, y, s = 2L, method = "mars", V = 4.0, number_of_bins = 10L)
#' regmdc(X_design, y, s = 2L, method = "mars", V = 4.0, number_of_bins = c(5L, 10L))
#' regmdc(X_design, y, s = 2L, method = "mars", V = 4.0, number_of_bins = c(NA, 10L))
#' regmdc(X_design, y, s = 2L, method = "tcmars", V = 2.0,
#'        constrained_interactions = c('1', '1-2'),
#'        increasing_interactions = c('2'))
#' @export
regmdc <- function(X_design, y, s, method, V = Inf, threshold = 1e-6,
                   number_of_bins = NULL,
                   constrained_interactions = NULL,
                   positive_interactions = NULL,
                   negative_interactions = NULL,
                   increasing_interactions = NULL,
                   decreasing_interactions = NULL) {
  # Error handling =============================================================
  if (length(dim(X_design)) != 2L) {
    stop('`X_design` must be a matrix or a data frame.')
  }
  if (anyNA(X_design)) {
    stop('`X_design` must not include NA or NaN.')
  }
  if (!all(apply(X_design, MARGIN = c(1L, 2L), FUN = is.numeric))) {
    stop('`X_design` must be numeric.')
  }
  if (min(X_design) < 0 || max(X_design) > 1) {
    stop('Every component of `X_design` must be between 0 and 1.')
  }

  if (anyNA(y)) {
    stop('`y` must not include NA or NaN.')
  }
  if (!is.numeric(y)) {
    stop('`y` must be numeric.')
  }
  if (length(y) != nrow(X_design)) {
    stop('`length(y)` must be equal to `nrow(X_design)`.')
  }

  if (!is.integer(s) || length(s) != 1L) {
    stop('`s` must be an integer.')
  }
  if (s < 1L || s > ncol(X_design)) {
    stop('`s` must be at least 1 and at most `ncol(X_design)`.')
  }

  if (!(method %in% c('em', 'hk', 'emhk', 'tc', 'mars', 'tcmars'))) {
    stop('`method` must be one of "em", "hk", "emhk", "tc", "mars", and "tcmars".')
  }

  if (!is.numeric(V) || length(V) != 1L) {
    stop('`V` must be a numeric scalar.')
  }
  if (V <= 0) {
    stop('`V` must be positive.')
  }

  if (!is.numeric(threshold) || length(threshold) != 1L) {
    stop('`threshold` must be a numeric scalar.')
  }
  if (threshold <= 0) {
    stop('`threshold` must be positive.')
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
      stop('Approximate methods are only available for `tc`, `mars`, and `tcmars` at this point.')
    }
  }
  # ============================================================================

  # Obtain the matrix for the LASSO problem and the vector indicating which
  # covariates are used in constructing each basis function. For totally convex
  # regression, MARS via LASSO, and their generalization, the indices of the
  # columns whose corresponding basis functions are constrained in the
  # estimation method are additionally collected.
  matrix_with_additional_info <- get_lasso_matrix(X_design, X_design, s,
                                                  method, number_of_bins)
  M <- matrix_with_additional_info$lasso_matrix
  basis_components <- matrix_with_additional_info$basis_components
  if (method %in% c('tc', 'mars', 'tcmars')) {
    constrained_basis <- matrix_with_additional_info$constrained_basis
  }
  is_included_basis <- matrix_with_additional_info$is_included_basis

  # Find the solution to the LASSO problem of the estimation method
  if (method == 'em') {
    constrained_cols <- NULL
    positive_cols <- (2L:ncol(M))
    negative_cols <- NULL
    solution <- solve_constrained_lasso(y, M, positive_indices = positive_cols)
  } else if (method == 'hk') {
    constrained_cols <- (2L:ncol(M))
    positive_cols <- NULL
    negative_cols <- NULL
    solution <- solve_constrained_lasso(y, M, V = V,
                                        constrained_indices = constrained_cols)
  } else if (method == 'emhk') {
    constrained_cols <- lapply(constrained_interactions, function(interaction) {
      which(basis_components == interaction)
    })
    positive_cols <- lapply(positive_interactions, function(interaction) {
      which(basis_components == interaction)
    })
    negative_cols <- lapply(negative_interactions, function(interaction) {
      which(basis_components == interaction)
    })
    constrained_cols <- unlist(constrained_cols)
    positive_cols <- unlist(positive_cols)
    negative_cols <- unlist(negative_cols)
    solution <- solve_constrained_lasso(y, M, V = V,
                                        constrained_indices = constrained_cols,
                                        positive_indices = positive_cols,
                                        negative_indices = negative_cols)
  } else if (method == 'tc') {
    constrained_cols <- NULL
    positive_cols <- constrained_basis
    negative_cols <- NULL
    solution <- solve_constrained_lasso(y, M, positive_indices = positive_cols)
  } else if (method == 'mars') {
    constrained_cols <- constrained_basis
    positive_cols <- NULL
    negative_cols <- NULL
    solution <- solve_constrained_lasso(y, M, V = V,
                                        constrained_indices = constrained_cols)
  } else {
    constrained_cols <- lapply(constrained_interactions, function(interaction) {
      utils::tail(which(basis_components == interaction), -1L)
    })
    positive_cols <- lapply(increasing_interactions, function(interaction) {
      utils::tail(which(basis_components == interaction), -1L)
    })
    negative_cols <- lapply(decreasing_interactions, function(interaction) {
      utils::tail(which(basis_components == interaction), -1L)
    })
    constrained_cols <- unlist(constrained_cols)
    positive_cols <- unlist(positive_cols)
    negative_cols <- unlist(negative_cols)
    solution <- solve_constrained_lasso(y, M, V = V,
                                        constrained_indices = constrained_cols,
                                        positive_indices = positive_cols,
                                        negative_indices = negative_cols)
  }

  names(solution) <- colnames(M)

  # Remove zero components and divide the solution into constrained and
  # unconstrained parts
  if (!is.null(constrained_cols)) {
    constrained_components <- solution[constrained_cols]
    is_nonzero_component <- (abs(constrained_components) >= threshold)
    constrained_components <- constrained_components[is_nonzero_component]
    V_solution <- sum(abs(constrained_components))

    unconstrained_components <- solution[-constrained_cols]
    is_nonzero_component <- (abs(unconstrained_components) >= threshold)
    unconstrained_components <- unconstrained_components[is_nonzero_component]
  } else {
    constrained_components <- NULL
    V_solution <- NULL
    unconstrained_components <- NULL
  }

  is_nonzero_component <- (abs(solution) >= threshold)
  compressed_solution <- solution[is_nonzero_component]

  # Compute the fitted values at the design points
  fitted_values <- compute_fit(X_design, X_design, s, method, number_of_bins,
                               compressed_solution, is_nonzero_component)

  # Compute the empirical loss (= mean squared error)
  empirical_loss <- compute_mse(y, fitted_values)
  # ============================================================================

  regmdc_model <- list(
    X_design = X_design,
    y = y,
    s = s,
    method = method,
    V = V,
    threshold = threshold,
    number_of_bins = number_of_bins,
    constrained_interactions = constrained_interactions,
    positive_interactions = positive_interactions,
    negative_interactions = negative_interactions,
    increasing_interactions = increasing_interactions,
    decreasing_interactions = decreasing_interactions,
    compressed_solution = compressed_solution,
    constrained_components = constrained_components,
    unconstrained_components = unconstrained_components,
    is_nonzero_component = is_nonzero_component,
    V_solution = V_solution,
    is_included_basis = is_included_basis,
    fitted_values = fitted_values,
    empirical_loss = empirical_loss
  )
  class(regmdc_model) <- "regmdc"

  regmdc_model
}
