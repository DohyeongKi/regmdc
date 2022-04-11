#' Construct the matrix for the LASSO problem of an estimation method
#'
#' Given an estimation method, this function constructs the matrix for the
#' corresponding LASSO problem. This function also returns the vector indicating
#' which covariates are used in constructing each basis function. Recall that
#' basis functions correspond to columns of the matrix for the problem. For
#' totally convex regression ("tc"), MARS via LASSO ("mars"), and their
#' generalization ("tcmars"), the indices of columns whose corresponding basis
#' functions are constrained in the estimation method are additionally returned.
#' See, for example, Section 3 of Fang et al. (2021) (for entirely monotonic
#' regression and Hardy—Krause variation denoising) and Section 3 and Section 6
#' of Ki et al. (2021) (for MARS via LASSO) for more details on the corresponding
#' LASSO problems. There are some other ongoing research and working papers, and
#' they would be available in the future.
#'
#' @param X_eval A numeric evaluation matrix. Each row corresponds to individual
#'   evaluation point at which basis functions are computed.
#' @param X_design A numeric design matrix. Each row corresponds to individual
#'   data. Basis functions are constructed from this matrix.
#' @param s A numeric scalar indicating the maximum order of interaction between
#'   covariates allowed in the estimation method.
#' @param method A string indicating the estimation method. One of "em", "hk",
#' "emhk", "tc", "mars", and "tcmars".
#' @references Ki, D., Fang, B., and Guntuboyina, A. (2021). MARS via LASSO.
#'   \url{https://arxiv.org/abs/2111.11694}.
#' @references Fang, B., Guntuboyina, A., and Sen, B. (2021). Multivariate
#' extensions of isotonic regression and total variation denoising via entire
#' monotonicity and Hardy—Krause variation. \emph{The Annals of Statistics},
#' \strong{49}(2), 769-792.
get_lasso_matrix <- function(X_eval, X_design, s, method) {
  if (method %in% c('em', 'hk', 'emhk')) {
    get_lasso_matrix_emhk(X_eval, X_design, s)
  } else if (method %in% c('tc', 'mars', 'tcmars')) {
    get_lasso_matrix_tcmars(X_eval, X_design, s)
  } else {
    stop('`method` must be one of "em", "hk", "emhk", "tc", "mars", and "tcmars".')
  }
}

get_lasso_matrix_emhk <- function(X_eval, X_design, s) {
  d <- ncol(X_design)
  unique_entries <- get_unique_column_entries(X_design, 'emhk')

  # Evaluate basis functions at the evaluation points ==========================
  lasso_matrix <- do.call(rbind, lapply((1L:nrow(X_eval)), function(row) {
    # Evaluate univariate indicator functions at the (row)-th evaluation point
    indicators <- lapply((1L:d), function(col) {
      c(TRUE, compute_indicator(X_eval[row, col] - unique_entries[[col]]))
    })

    # Record the order of each univariate indicator function. We say constant
    # functions have order zero and univariate indicator functions have order one.
    indicators_orders <- lapply((1L:d), function(col) {
      c(0L, rep(1L, length(unique_entries[[col]])))
    })

    # Evaluate basis functions at the (row)-th evaluation point
    basis <- indicators[[1L]]
    # Compute the order of each basis function which is defined as the number of
    # univariate indicator functions multiplied
    basis_orders <- indicators_orders[[1L]]
    # Put aside basis functions that reach the maximum allowed order s in the
    # for loop below to reduce computational cost
    full_order_basis <- c()

    for (k in (2L:d)) {
      basis <- outer(basis, indicators[[k]], "&")
      basis <- c(basis)
      basis_orders <- outer(basis_orders, indicators_orders[[k]], "+")
      basis_orders <- c(basis_orders)

      is_full_order <- (basis_orders == s)
      full_order_basis <- append(full_order_basis, basis[is_full_order])
      basis <- basis[!is_full_order]
      basis_orders <- basis_orders[!is_full_order]
    }

    basis <- append(basis, full_order_basis)
  }))

  lasso_matrix <- apply(lasso_matrix, MARGIN = c(1, 2), FUN = as.numeric)
  # ============================================================================

  # Create column names of the matrix ==========================================
  # Give each univariate indicator function a name
  indicators_names <- lapply((1L:d), function(col) {
    names <- sapply((1L:length(unique_entries[[col]])), function(k) {
      paste0(paste0("I(Var", col, "-",
                    format(unique_entries[[col]][k], digits = 4L), ")"),
             collapse = "*")
    })
    c("", names)
  })

  # Record the order of each univariate indicator function. We say constant
  # functions have order zero and univariate indicator functions have order one.
  indicators_orders <- lapply((1L:d), function(col) {
    c(0L, rep(1L, length(unique_entries[[col]])))
  })

  # Record the covariate from which each univariate indicator function is
  # constructed
  indicators_components <- lapply((1L:d), function(col) {
    c("", rep(toString(col), length(unique_entries[[col]])))
  })

  # Give each basis function a name
  basis_names <- indicators_names[[1L]]
  # Compute the order of each basis function which is defined as the number of
  # univariate indicator functions multiplied
  basis_orders <- indicators_orders[[1L]]
  # Find the covariates from which each basis function is constructed
  basis_components <- indicators_components[[1L]]

  # Put aside basis functions that reach the maximum allowed order s in the
  # for loop to reduce computational cost
  full_order_basis_names <- c()
  full_order_basis_components <- c()

  for (k in (2L:d)) {
    basis_names <- outer(basis_names, indicators_names[[k]], "paste0")
    basis_names <- c(basis_names)
    basis_orders <- outer(basis_orders, indicators_orders[[k]], "+")
    basis_orders <- c(basis_orders)
    basis_components <- outer(basis_components, indicators_components[[k]],
                              Vectorize(paste_with_hyphen))
    basis_components <- c(basis_components)

    is_full_order <- (basis_orders == s)
    full_order_basis_names <- append(full_order_basis_names,
                                     basis_names[is_full_order])
    full_order_basis_components <- append(full_order_basis_components,
                                          basis_components[is_full_order])
    basis_names <- basis_names[!is_full_order]
    basis_orders <- basis_orders[!is_full_order]
    basis_components <- basis_components[!is_full_order]
  }

  basis_names <- append(basis_names, full_order_basis_names)
  basis_names[1L] <- "(Intercept)"
  colnames(lasso_matrix) <- basis_names

  basis_components <- append(basis_components, full_order_basis_components)
  # ============================================================================

  list(
    lasso_matrix = lasso_matrix,
    basis_components = basis_components
  )
}


get_lasso_matrix_tcmars <- function(X_eval, X_design, s) {
  d <- ncol(X_design)
  unique_entries <- get_unique_column_entries(X_design, 'tcmars')

  # Evaluate basis functions at the evaluation points ==========================
  lasso_matrix <- do.call(rbind, lapply((1L:nrow(X_eval)), function(row) {
    # Evaluate hinge functions at the (row)-th evaluation point
    hinges <- lapply((1L:d), function(col) {
      c(1L, compute_hinge(X_eval[row, col] - unique_entries[[col]]))
    })

    # Record the order of each hinge function. We say constant functions have
    # order zero and hinge functions have order one.
    hinges_orders <- lapply((1L:d), function(col) {
      c(0L, rep(1L, length(unique_entries[[col]])))
    })

    # Evaluate basis functions at the (row)-th evaluation point
    basis <- hinges[[1L]]
    # Compute the order of each basis function which is defined as the number of
    # hinge functions multiplied
    basis_orders <- hinges_orders[[1L]]
    # Put aside basis functions that reach the maximum allowed order s in the
    # for loop below to reduce computational cost
    full_order_basis <- c()

    for (k in (2L:d)) {
      basis <- outer(basis, hinges[[k]], "*")
      basis <- c(basis)
      basis_orders <- outer(basis_orders, hinges_orders[[k]], "+")
      basis_orders <- c(basis_orders)

      is_full_order <- (basis_orders == s)
      full_order_basis <- append(full_order_basis, basis[is_full_order])
      basis <- basis[!is_full_order]
      basis_orders <- basis_orders[!is_full_order]
    }

    basis <- append(basis, full_order_basis)
  }))
  # ============================================================================

  # Create column names and find the constrained columns of the matrix =========
  # Give each hinge function a name
  hinges_names <- lapply((1L:d), function(col) {
    names <- sapply((1L:length(unique_entries[[col]])), function(k) {
      paste0(paste0("H(Var", col, "-",
                    format(unique_entries[[col]][k], digits = 4L), ")"),
             collapse = "*")
    })
    c("", names)
  })

  # Record the order of each hinge function. We say constant functions have
  # order zero and hinge functions have order one.
  hinges_orders <- lapply((1L:d), function(col) {
    c(0L, rep(1L, length(unique_entries[[col]])))
  })

  # Record the covariate from which each hinge function is constructed
  hinges_components <- lapply((1L:d), function(col) {
    c("", rep(toString(col), length(unique_entries[[col]])))
  })

  # Record whether each hinge function is constrained or not. We do not
  # constrain constant functions and linear functions.
  is_constrained_hinge <- lapply((1L:d), function(col) {
    c(FALSE, FALSE, rep(TRUE, (length(unique_entries[[col]]) - 1L)))
  })

  # Give each basis function a name
  basis_names <- hinges_names[[1L]]
  # Compute the order of each basis function which is defined as the number of
  # hinge functions multiplied
  basis_orders <- hinges_orders[[1L]]
  # Find the covariates from which each basis function is constructed
  basis_components <- hinges_components[[1L]]
  # Check whether each basis function is constrained or not
  is_constrained_basis <- is_constrained_hinge[[1L]]

  # Put aside basis functions that reach the maximum allowed order s in the
  # for loop below to reduce computational cost
  full_order_basis_names <- c()
  full_order_basis_components <- c()
  is_constrained_full_order_basis <- c()

  for (k in (2L:d)) {
    basis_names <- outer(basis_names, hinges_names[[k]], "paste0")
    basis_names <- c(basis_names)
    basis_orders <- outer(basis_orders, hinges_orders[[k]], "+")
    basis_orders <- c(basis_orders)
    basis_components <- outer(basis_components, hinges_components[[k]],
                              Vectorize(paste_with_hyphen))
    basis_components <- c(basis_components)
    is_constrained_basis <- outer(is_constrained_basis,
                                  is_constrained_hinge[[k]], "|")
    is_constrained_basis <- c(is_constrained_basis)

    is_full_order <- (basis_orders == s)
    full_order_basis_names <- append(full_order_basis_names,
                                     basis_names[is_full_order])
    full_order_basis_components <- append(full_order_basis_components,
                                          basis_components[is_full_order])
    is_constrained_full_order_basis <- append(is_constrained_full_order_basis,
                                              is_constrained_basis[is_full_order])
    basis_names <- basis_names[!is_full_order]
    basis_orders <- basis_orders[!is_full_order]
    basis_components <- basis_components[!is_full_order]
    is_constrained_basis <- is_constrained_basis[!is_full_order]
  }

  basis_names <- append(basis_names, full_order_basis_names)
  basis_names[1L] <- "(Intercept)"
  colnames(lasso_matrix) <- basis_names

  basis_components <- append(basis_components, full_order_basis_components)

  is_constrained_basis <- append(is_constrained_basis,
                                 is_constrained_full_order_basis)
  constrained_basis <- which(is_constrained_basis == TRUE)
  # ============================================================================

  list(
    lasso_matrix = lasso_matrix,
    basis_components = basis_components,
    constrained_basis = constrained_basis
  )
}
