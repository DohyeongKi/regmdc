#' Construct the matrix for the LASSO problem of an estimation method
#'
#' Given an estimation method, this function constructs the matrix for the
#' corresponding LASSO problem. This function also returns the vector indicating
#' which covariates are used in constructing each basis function. Recall that
#' basis functions correspond to columns of the matrix for the problem. For
#' totally convex regression ("tc"), MARS via LASSO ("mars"), and their
#' generalization ("tcmars"), the indices of columns whose corresponding basis
#' functions are constrained in the estimation method are additionally returned.
#'
#' @param X_eval A numeric evaluation matrix. Each row corresponds to an
#'   individual evaluation point at which basis functions are computed.
#' @param X_design A numeric design matrix. Each row corresponds to an
#'   individual data. Basis functions are constructed from this matrix.
#' @param s A numeric scalar indicating the maximum order of interaction between
#'   covariates allowed in the estimation method.
#' @param method A string indicating the estimation method. One of "em", "hk",
#'   "emhk", "tc", "mars", and "tcmars".
#' @param is_lattice A logical scalar for whether the design is lattice or not.
#'   Only used for "em", "hk", and "emhk".
#' @param number_of_bins An integer vector of the numbers of bins for the
#'   approximate methods. `NULL` if the approximate methods are not used.
#'   Currently available for "tc", "mars", and "tcmars".
#' @param is_included_basis A logical vector indicating whether or not each
#'   basis function is included in the LASSO problem.
#' @references Ki, D., Fang, B., and Guntuboyina, A. (2021). MARS via LASSO.
#'   \url{https://arxiv.org/abs/2111.11694}.
#' @references Fang, B., Guntuboyina, A., and Sen, B. (2021). Multivariate
#'   extensions of isotonic regression and total variation denoising via entire
#'   monotonicity and Hardy—Krause variation. \emph{The Annals of Statistics},
#'   \strong{49}(2), 769-792.
get_lasso_matrix <- function(X_eval, X_design, s, method, is_lattice,
                             number_of_bins, is_included_basis = NULL) {
  if (method %in% c('em', 'hk', 'emhk')) {
    if (is_lattice) {
      get_lasso_matrix_emhk_lattice(X_eval, X_design, s, is_included_basis)
    } else {
      get_lasso_matrix_emhk_nonlattice(X_eval, X_design, s, is_included_basis)
    }
  } else if (method %in% c('tc', 'mars', 'tcmars')) {
    get_lasso_matrix_tcmars(X_eval, X_design, s, number_of_bins,
                            is_included_basis)
  } else {
    stop('`method` must be one of "em", "hk", "emhk", "tc", "mars", and "tcmars".')
  }
}


get_lasso_matrix_emhk_lattice <- function(X_eval, X_design, s,
                                          is_included_basis = NULL) {
  d <- ncol(X_design)
  unique_entries <- get_unique_column_entries(X_design, 'emhk')
  for (col in (1L:d)) {
    if (length(unique_entries[[col]]) == 0L) {
      if (is.null(colnames(X_design))) {
        stop(paste0('All the values of Var', col, ' are zero. Please remove that variable.'))
      } else {
        stop(paste0('All the values of "', colnames(X_design)[col], '" are zero. Please remove that variable.'))
      }
    }
  }

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

    if (d >= 2L) {
      for (k in (2L:d)) {
        is_full_order <- (basis_orders == s)
        full_order_basis <- append(full_order_basis, basis[is_full_order])
        basis <- basis[!is_full_order]
        basis_orders <- basis_orders[!is_full_order]

        basis <- outer(basis, indicators[[k]], "&")
        basis <- c(basis)
        basis_orders <- outer(basis_orders, indicators_orders[[k]], "+")
        basis_orders <- c(basis_orders)
      }
    }

    basis <- append(basis, full_order_basis)
  }))

  lasso_matrix <- apply(lasso_matrix, MARGIN = c(1, 2), FUN = as.numeric)
  # ============================================================================

  # Create column names of the matrix ==========================================
  # Give a name to each univariate indicator function
  if (is.null(colnames(X_design))) {
    indicators_names <- lapply((1L:d), function(col) {
      names <- sapply((1L:length(unique_entries[[col]])), function(k) {
        paste0("I(Var", col, "-",
               format(unique_entries[[col]][k], digits = 4L), ")")
      })
      c("", names)
    })
  } else {
    indicators_names <- lapply((1L:d), function(col) {
      names <- sapply((1L:length(unique_entries[[col]])), function(k) {
        paste0("I(", colnames(X_design)[col], "-",
               format(unique_entries[[col]][k], digits = 4L), ")")
      })
      c("", names)
    })
  }

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

  # Give a name to each basis function
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

  if (d >= 2L) {
    for (k in (2L:d)) {
      is_full_order <- (basis_orders == s)
      full_order_basis_names <- append(full_order_basis_names,
                                       basis_names[is_full_order])
      full_order_basis_components <- append(full_order_basis_components,
                                            basis_components[is_full_order])
      basis_names <- basis_names[!is_full_order]
      basis_orders <- basis_orders[!is_full_order]
      basis_components <- basis_components[!is_full_order]

      basis_names <- outer(basis_names, indicators_names[[k]], "paste0")
      basis_names <- c(basis_names)
      basis_orders <- outer(basis_orders, indicators_orders[[k]], "+")
      basis_orders <- c(basis_orders)
      basis_components <- outer(basis_components, indicators_components[[k]],
                                Vectorize(paste_with_hyphen))
      basis_components <- c(basis_components)
    }

    basis_names <- append(basis_names, full_order_basis_names)
    basis_components <- append(basis_components, full_order_basis_components)
  }

  basis_names[1L] <- "(Intercept)"
  colnames(lasso_matrix) <- basis_names

  if (is.null(is_included_basis)) {
    # Remove all zero columns
    is_included_basis <- apply((lasso_matrix != 0), MARGIN = 2, any)
  } else {
    if (length(is_included_basis) != ncol(lasso_matrix)) {
      stop('`length(is_included_basis)` must be equal to the number of columns of the LASSO matrix.')
    }
  }

  # Include specified columns
  lasso_matrix <- lasso_matrix[, is_included_basis]
  basis_components <- basis_components[is_included_basis]
  # ============================================================================

  list(
    lasso_matrix = lasso_matrix,
    basis_components = basis_components,
    is_included_basis = is_included_basis
  )
}


get_lasso_matrix_emhk_nonlattice <- function(X_eval, X_design, s,
                                             is_included_basis = NULL) {
  d <- ncol(X_design)
  unique_entries <- get_unique_column_entries(X_design, 'emhk')
  for (col in (1L:d)) {
    if (length(unique_entries[[col]]) == 0L) {
      if (is.null(colnames(X_design))) {
        stop(paste0('All the values of Var', col, ' are zero. Please remove that variable.'))
      } else {
        stop(paste0('All the values of "', colnames(X_design)[col], '" are zero. Please remove that variable.'))
      }
    }
  }

  # Consider the constant term =================================================
  lasso_matrix <- matrix(rep(1.0, nrow(X_eval)), ncol = 1L)
  # Give a name to the basis function
  colnames(lasso_matrix) <- c("(Intercept)")
  # Record the covariates from which the function is constructed
  basis_components <- c("")
  # ============================================================================

  # Add the first order terms ==================================================
  for (col in (1L:d)) {
    column_unique <- unique_entries[[col]]

    X_eval_col <- X_eval[, col]
    basis <- sapply(column_unique, simplify = TRUE, function(entry) {
      compute_indicator(X_eval_col - entry)
    })

    if (is.null(colnames(X_design))) {
      basis_names <- sapply(column_unique, simplify = TRUE, function(entry) {
        paste0("I(Var", col, "-", format(entry, digits = 4L), ")")
      })
      colnames(basis) <- basis_names
    } else {
      basis_names <- sapply(column_unique, simplify = TRUE, function(entry) {
        paste0("I(", colnames(X_design)[col], "-", format(entry, digits = 4L), ")")
      })
      colnames(basis) <- basis_names
    }

    lasso_matrix <- cbind(lasso_matrix, basis)
    basis_components <- c(basis_components,
                          rep(toString(col), length(column_unique)))
  }
  # ============================================================================

  # Add the higher order terms =================================================
  if (s >= 2L) {
    for (order in (2L:s)) {
      subsets <- utils::combn(d, s)  # all subsets of {1, ... , d} of size s

      for (index in ncol(subsets)) {
        subset <- subsets[, index]  # a particular subset of size s
        # Only consider the corresponding columns
        X_subset <- X_design[, subset]
        nonzero_rows <- apply(X_subset, MARGIN = 1L, function(row) {
          all(row > 0)
        })
        X_subset <- X_subset[nonzero_rows, ]
        X_subset <- rbind(X_subset, rep(1.0, s))

        # Compute component-wise minimum vectors
        X_min <- X_subset
        for (iter in (1L:s)) {
          X_min_new <- do.call(rbind, lapply((1L:nrow(X_subset)), function(row) {
            sapply((1L:s), simplify = TRUE, function(col) {
              pmin(X_min[, col], X_subset[row, col])
            })
          }))
          colnames(X_min_new) <- colnames(X_min)

          X_min <- rbind(X_min, X_min_new)
          X_min <- unique(X_min, MARGIN = 1L)
        }
        X_min <- utils::head(X_min, -1)

        # Evaluate the basis functions at the evaluation points
        X_eval_subset <- X_eval[, subset]
        basis <- sapply((1L:nrow(X_min)), simplify = TRUE, function(row) {
          X_min_row <- X_min[row, ]
          apply(X_eval_subset, MARGIN = 1L, FUN = function(eval_point) {
            all(compute_indicator(eval_point - X_min_row))
          })
        })

        # Give names to the basis functions
        if (is.null(colnames(X_design))) {
          basis_names <- sapply((1L:nrow(X_min)), simplify = TRUE, function(row) {
            basis_name <- ""
            for (col in (1L:ncol(X_min))) {
              basis_name <- paste0(basis_name,
                                   "I(Var", col, "-",
                                   format(X_min[row, col], digits = 4L), ")"
              )
            }
            basis_name
          })
          colnames(basis) <- basis_names
        } else {
          basis_names <- sapply((1L:nrow(X_min)), simplify = TRUE, function(row) {
            basis_name <- ""
            for (col in (1L:ncol(X_min))) {
              basis_name <- paste0(basis_name,
                                   "I(", colnames(X_design)[col], "-",
                                   format(X_min[row, col], digits = 4L), ")")
            }
            basis_name
          })
          colnames(basis) <- basis_names
        }

        # Find the covariates from which the basis functions are constructed
        basis_component <- subset[1]
        if (length(subset) >= 2L) {
          for (iter in (2L:length(subset))) {
            basis_component <- paste(basis_component, subset[iter], sep = "-")
          }
        }

        lasso_matrix <- cbind(lasso_matrix, basis)
        basis_components <- c(basis_components, rep(basis_component, nrow(X_min)))
      }
    }
  }

  if (is.null(is_included_basis)) {
    is_included_basis <- rep(TRUE, ncol(lasso_matrix))
  } else {
    if (length(is_included_basis) != ncol(lasso_matrix)) {
      stop('`length(is_included_basis)` must be equal to the number of columns of the LASSO matrix.')
    }
  }

  # Include specified columns
  lasso_matrix <- lasso_matrix[, is_included_basis]
  basis_components <- basis_components[is_included_basis]
  # ============================================================================

  list(
    lasso_matrix = lasso_matrix,
    basis_components = basis_components,
    is_included_basis = is_included_basis
  )
  # ============================================================================
}


get_lasso_matrix_tcmars <- function(X_eval, X_design, s, number_of_bins,
                                    is_included_basis = NULL) {
  d <- ncol(X_design)
  unique_entries <- get_unique_column_entries(X_design, 'tcmars', number_of_bins)
  for (col in (1L:d)) {
    if (length(unique_entries[[col]]) == 0L) {
      if (is.null(colnames(X_design))) {
        stop(paste0('All the values of Var', col, ' are zero. Please remove that variable.'))
      } else {
        stop(paste0('All the values of "', colnames(X_design)[col], '" are zero. Please remove that variable.'))
      }
    }
  }

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

    if (d >= 2L) {
      for (k in (2L:d)) {
        is_full_order <- (basis_orders == s)
        full_order_basis <- append(full_order_basis, basis[is_full_order])
        basis <- basis[!is_full_order]
        basis_orders <- basis_orders[!is_full_order]

        basis <- outer(basis, hinges[[k]], "*")
        basis <- c(basis)
        basis_orders <- outer(basis_orders, hinges_orders[[k]], "+")
        basis_orders <- c(basis_orders)
      }
    }

    basis <- append(basis, full_order_basis)
  }))
  # ============================================================================

  # Create column names and find the constrained columns of the matrix =========
  # Give a name to each hinge function
  if (is.null(colnames(X_design))) {
    hinges_names <- lapply((1L:d), function(col) {
      names <- sapply((1L:length(unique_entries[[col]])), function(k) {
        paste0("H(Var", col, "-",
               format(unique_entries[[col]][k], digits = 4L), ")")
      })
      c("", names)
    })
  } else {
    hinges_names <- lapply((1L:d), function(col) {
      names <- sapply((1L:length(unique_entries[[col]])), function(k) {
        paste0("H(", colnames(X_design)[col], "-",
               format(unique_entries[[col]][k], digits = 4L), ")")
      })
      c("", names)
    })
  }

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

  # Give a name to each basis function
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

  if (d >= 2L) {
    for (k in (2L:d)) {
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
    }

    basis_names <- append(basis_names, full_order_basis_names)
    basis_components <- append(basis_components, full_order_basis_components)
    is_constrained_basis <- append(is_constrained_basis,
                                   is_constrained_full_order_basis)
  }

  basis_names[1L] <- "(Intercept)"
  colnames(lasso_matrix) <- basis_names

  if (is.null(is_included_basis)) {
    # Remove all zero columns
    is_included_basis <- apply((lasso_matrix != 0), MARGIN = 2, any)
  } else {
    if (length(is_included_basis) != ncol(lasso_matrix)) {
      stop('`length(is_included_basis)` must be equal to the number of columns of the LASSO matrix.')
    }
  }

  # Include specified columns
  lasso_matrix <- lasso_matrix[, is_included_basis]
  basis_components <- basis_components[is_included_basis]
  is_constrained_basis <- is_constrained_basis[is_included_basis]
  constrained_basis <- which(is_constrained_basis == TRUE)
  # ============================================================================

  list(
    lasso_matrix = lasso_matrix,
    basis_components = basis_components,
    constrained_basis = constrained_basis,
    is_included_basis = is_included_basis
  )
}
