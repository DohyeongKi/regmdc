#' Construct the matrix for the LASSO problem of an estimation method
#'
#' Given an estimation method, this function constructs the matrix for the
#' corresponding LASSO problem. This function also returns the logical vectors
#' for whether each basis function is restricted to have a positive coefficient
#' or a negative coefficient and whether its variation is constrained or not in
#' the model. Recall that basis functions correspond to columns of the matrix
#' for the LASSO problem. For totally concave regression ("tc"), MARS via LASSO
#' ("mars"), and their generalization ("tcmars"), the scale factors of basis
#' functions, which are needed for rescaling, are additionally returned.
#'
#' @param X_eval A numeric evaluation matrix. Each row corresponds to an
#'   individual evaluation point at which basis functions are computed.
#' @param X_design A numeric design matrix. Each row corresponds to an
#'   individual data. Basis functions are constructed from this matrix.
#' @param s A numeric scalar indicating the maximum order of interaction between
#'   covariates allowed in the estimation method.
#' @param method A string indicating the estimation method. One of "em", "hk",
#'   "emhk", "tc", "mars", and "tcmars".
#' @param is_scaled A logical scalar for whether the design matrix is scaled so
#'   that every entry is between 0 and 1. If `FALSE`, the min-max scaling is
#'   applied to each column of the design matrix.
#' @param is_lattice A logical scalar for whether the design is lattice or not.
#'   Only used for "em", "hk", and "emhk".
#' @param number_of_bins An integer or an integer vector of the numbers of bins
#'   for the approximate method. Currently available for "tc", "mars", and
#'   "tcmars".
#' @param increasing_covariates An integer vector of monotonically increasing
#'   covariates. Possibly used for "em" and "emhk".
#' @param decreasing_covariates An integer vector of monotonically decreasing
#'   covariates. Possibly used for "em" and "emhk".
#' @param concave_covariates An integer vector of concave covariates. Possibly
#'   used for "tc" and "tcmars".
#' @param convex_covariates An integer vector of convex covariates. Possibly
#'   used for "tc" and "tcmars".
#' @param variation_constrained_covariates An integer vector of covariates whose
#'   variation is constrained. Used for "hk" and "mars" and also possibly used
#'   for "emhk" and "tcmars".
#' @param extra_linear_covariates An integer vector or a string vector of extra
#'   linear covariates added to the model. Possibly used for "tc", "mars", and
#'   "tcmars".
#' @param is_included_basis A logical vector indicating whether or not each
#'   basis function is included in the LASSO problem.
#' @references Ki, D., Fang, B., and Guntuboyina, A. (2024+). MARS via LASSO.
#'   Accepted at \emph{Annals of Statistics}. Available at
#'   \url{https://arxiv.org/abs/2111.11694}.
#' @references Fang, B., Guntuboyina, A., and Sen, B. (2021). Multivariate
#'   extensions of isotonic regression and total variation denoising via entire
#'   monotonicity and Hardy—Krause variation. \emph{Annals of Statistics},
#'   \strong{49}(2), 769-792.
get_lasso_matrix <- function(X_eval, X_design, s, method, is_scaled, is_lattice,
                             number_of_bins,
                             increasing_covariates, decreasing_covariates,
                             concave_covariates, convex_covariates,
                             variation_constrained_covariates,
                             extra_linear_covariates, is_included_basis = NULL) {
  # Give names to the columns of the design matrix if there aren't
  if (is.null(colnames(X_design))) {
    colnames(X_design) <- paste0("Var", (1L:ncol(X_design)))
  }

  # Scale the matrices if necessary. Record the maximal and minimal values of
  # each column for scaling back entries later.
  if (is_scaled) {
    max_vals <- rep(1, ncol(X_design))
    min_vals <- rep(0, ncol(X_design))
  } else {
    max_vals <- apply(X_design, MARGIN = 2L, max)
    min_vals <- apply(X_design, MARGIN = 2L, min)

    for (col in (1L:ncol(X_design))) {
      if (max_vals[col] == min_vals[col]) {
        stop(paste0('All the values of "', colnames(X_design)[col], '" are the same. Please remove that variable.'))
      } else {
        X_design[, col] <- ((X_design[, col] - min_vals[col])
                            / (max_vals[col] - min_vals[col]))
        X_eval[, col] <- ((X_eval[, col] - min_vals[col])
                          / (max_vals[col] - min_vals[col]))
      }
    }
  }

  if (method %in% c('em', 'hk', 'emhk')) {
    if (is_lattice) {
      get_lasso_matrix_emhk_lattice(
        X_eval, X_design, max_vals, min_vals, s,
        increasing_covariates, decreasing_covariates,
        variation_constrained_covariates, is_included_basis
      )
    } else {
      get_lasso_matrix_emhk_nonlattice(
        X_eval, X_design, max_vals, min_vals, s,
        increasing_covariates, decreasing_covariates,
        variation_constrained_covariates, is_included_basis
      )
    }
  } else if (method %in% c('tc', 'mars', 'tcmars')) {
    get_lasso_matrix_tcmars(
      X_eval, X_design, max_vals, min_vals, s, number_of_bins,
      concave_covariates, convex_covariates,
      variation_constrained_covariates, extra_linear_covariates,
      is_included_basis
    )
  } else {
    stop('`method` must be one of "em", "hk", "emhk", "tc", "mars", and "tcmars".')
  }
}


get_lasso_matrix_emhk_lattice <- function(X_eval, X_design,
                                          max_vals, min_vals, s,
                                          increasing_covariates,
                                          decreasing_covariates,
                                          variation_constrained_covariates,
                                          is_included_basis = NULL) {
  # Record whether basis dropping can happen due to the sign constraints
  is_basis_drop_possible <- !(is.null(increasing_covariates)
                              || is.null(decreasing_covariates))

  d <- ncol(X_design)
  unique_entries <- get_unique_column_entries(X_design, 'emhk')
  for (col in (1L:d)) {
    if (length(unique_entries[[col]]) == 0L) {
      stop(paste0('All the values of "', colnames(X_design)[col], '" are zero. Please remove that variable.'))
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

    if (is_basis_drop_possible) {
      # Record whether each univariate indicator function is restricted to have
      # a positive coefficient or a negative coefficient and whether its
      # variation is constrained or not
      is_positive_indicator <- lapply((1L:d), function(col) {
        if (col %in% increasing_covariates) {
          c(FALSE, rep(TRUE, length(unique_entries[[col]])))
        } else {
          rep(FALSE, length(unique_entries[[col]]) + 1L)
        }
      })
      is_negative_indicator <- lapply((1L:d), function(col) {
        if (col %in% decreasing_covariates) {
          c(FALSE, rep(TRUE, length(unique_entries[[col]])))
        } else {
          rep(FALSE, length(unique_entries[[col]]) + 1L)
        }
      })
      is_variation_constrained_indicator <- lapply((1L:d), function(col) {
        if (col %in% variation_constrained_covariates) {
          c(FALSE, rep(TRUE, length(unique_entries[[col]])))
        } else {
          rep(FALSE, length(unique_entries[[col]]) + 1L)
        }
      })
    }

    # Evaluate basis functions at the (row)-th evaluation point
    basis <- indicators[[1L]]
    # Compute the order of each basis function which is defined as the number of
    # univariate indicator functions multiplied
    basis_orders <- indicators_orders[[1L]]
    if (is_basis_drop_possible) {
      # Record whether each basis function is restricted to have a positive
      # coefficient or a negative coefficient and whether its variation is
      # constrained or not
      is_positive_basis <- is_positive_indicator[[1L]]
      is_negative_basis <- is_negative_indicator[[1L]]
      is_variation_constrained_basis <- is_variation_constrained_indicator[[1L]]
    }
    # Put aside basis functions that reach the maximum allowed order s in the
    # for loop below to reduce computational cost
    full_order_basis <- c()

    if (d >= 2L) {
      for (k in (2L:d)) {
        # Put aside basis functions that reach the maximum allowed order s
        is_full_order <- (basis_orders == s)
        full_order_basis <- append(full_order_basis, basis[is_full_order])
        basis <- basis[!is_full_order]
        basis_orders <- basis_orders[!is_full_order]

        basis <- outer(basis, indicators[[k]], "&")
        basis <- c(basis)
        basis_orders <- outer(basis_orders, indicators_orders[[k]], "+")
        basis_orders <- c(basis_orders)

        if (is_basis_drop_possible) {
          is_positive_basis <- is_positive_basis[!is_full_order]
          is_negative_basis <- is_negative_basis[!is_full_order]
          is_variation_constrained_basis <- (
            is_variation_constrained_basis[!is_full_order]
          )

          # The descendants of basis functions that have positive (resp.
          # negative) coefficients are required to have positive (resp. negative)
          # coefficients. The descendants of basis functions whose variation is
          # constrained are also under the variation constraint.
          is_positive_basis <- outer(is_positive_basis,
                                     is_positive_indicator[[k]], "|")
          is_positive_basis <- c(is_positive_basis)
          is_negative_basis <- outer(is_negative_basis,
                                     is_negative_indicator[[k]], "|")
          is_negative_basis <- c(is_negative_basis)
          is_variation_constrained_basis <- outer(
            is_variation_constrained_basis,
            is_variation_constrained_indicator[[k]], "|"
          )
          is_variation_constrained_basis <- c(is_variation_constrained_basis)

          # Remove basis functions whose coefficients need to be both positive
          # and negative
          is_nonzero <- !(is_positive_basis & is_negative_basis)
          basis <- basis[is_nonzero]
          basis_orders <- basis_orders[is_nonzero]
          is_positive_basis <- is_positive_basis[is_nonzero]
          is_negative_basis <- is_negative_basis[is_nonzero]
          is_variation_constrained_basis <- (
            is_variation_constrained_basis[is_nonzero]
          )
        }
      }
    }

    basis <- append(basis, full_order_basis)
  }))

  lasso_matrix <- apply(lasso_matrix, MARGIN = c(1, 2), FUN = as.numeric)
  # ============================================================================

  # Create column names of the matrix ==========================================
  # Give a name to each univariate indicator function
  indicators_names <- lapply((1L:d), function(col) {
    names <- sapply((1L:length(unique_entries[[col]])), function(k) {
      paste0("I(", colnames(X_design)[col], "-",
             scale_back_matrix_entry(unique_entries[[col]][k],
                                     max_vals[col], min_vals[col],
                                     digits = 4L),
             ")")
    })
    c("", names)
  })

  # Record the order of each univariate indicator function. We say constant
  # functions have order zero and univariate indicator functions have order one.
  indicators_orders <- lapply((1L:d), function(col) {
    c(0L, rep(1L, length(unique_entries[[col]])))
  })

  # Record whether each univariate indicator function is restricted to have a
  # positive coefficient or a negative coefficient and whether its variation
  # is constrained or not
  is_positive_indicator <- lapply((1L:d), function(col) {
    if (col %in% increasing_covariates) {
      c(FALSE, rep(TRUE, length(unique_entries[[col]])))
    } else {
      rep(FALSE, length(unique_entries[[col]]) + 1L)
    }
  })
  is_negative_indicator <- lapply((1L:d), function(col) {
    if (col %in% decreasing_covariates) {
      c(FALSE, rep(TRUE, length(unique_entries[[col]])))
    } else {
      rep(FALSE, length(unique_entries[[col]]) + 1L)
    }
  })
  is_variation_constrained_indicator <- lapply((1L:d), function(col) {
    if (col %in% variation_constrained_covariates) {
      c(FALSE, rep(TRUE, length(unique_entries[[col]])))
    } else {
      rep(FALSE, length(unique_entries[[col]]) + 1L)
    }
  })

  # Give a name to each basis function
  basis_names <- indicators_names[[1L]]
  # Compute the order of each basis function which is defined as the number of
  # univariate indicator functions multiplied
  basis_orders <- indicators_orders[[1L]]
  # Record whether each basis function is restricted to have a positive
  # coefficient or a negative coefficient and whether its variation is
  # constrained or not
  is_positive_basis <- is_positive_indicator[[1L]]
  is_negative_basis <- is_negative_indicator[[1L]]
  is_variation_constrained_basis <- is_variation_constrained_indicator[[1L]]

  # Put aside basis functions that reach the maximum allowed order s in the
  # for loop to reduce computational cost
  full_order_basis_names <- c()
  is_positive_full_order_basis <- c()
  is_negative_full_order_basis <- c()
  is_variation_constrained_full_order_basis <- c()

  if (d >= 2L) {
    for (k in (2L:d)) {
      # Put aside basis functions that reach the maximum allowed order s
      is_full_order <- (basis_orders == s)
      full_order_basis_names <- append(full_order_basis_names,
                                       basis_names[is_full_order])
      is_positive_full_order_basis <- append(
        is_positive_full_order_basis, is_positive_basis[is_full_order]
      )
      is_negative_full_order_basis <- append(
        is_negative_full_order_basis, is_negative_basis[is_full_order]
      )
      is_variation_constrained_full_order_basis <- append(
        is_variation_constrained_full_order_basis,
        is_variation_constrained_basis[is_full_order]
      )
      basis_names <- basis_names[!is_full_order]
      basis_orders <- basis_orders[!is_full_order]
      is_positive_basis <- is_positive_basis[!is_full_order]
      is_negative_basis <- is_negative_basis[!is_full_order]
      is_variation_constrained_basis <- (
        is_variation_constrained_basis[!is_full_order]
      )

      basis_names <- outer(basis_names, indicators_names[[k]], "paste0")
      basis_names <- c(basis_names)
      basis_orders <- outer(basis_orders, indicators_orders[[k]], "+")
      basis_orders <- c(basis_orders)

      # The descendants of basis functions that have positive (resp.
      # negative) coefficients are required to have positive (resp. negative)
      # coefficients. The descendants of basis functions whose variation is
      # constrained are also under the variation constraint
      is_positive_basis <- outer(is_positive_basis,
                                 is_positive_indicator[[k]], "|")
      is_positive_basis <- c(is_positive_basis)
      is_negative_basis <- outer(is_negative_basis,
                                 is_negative_indicator[[k]], "|")
      is_negative_basis <- c(is_negative_basis)
      is_variation_constrained_basis <- outer(
        is_variation_constrained_basis,
        is_variation_constrained_indicator[[k]], "|"
      )
      is_variation_constrained_basis <- c(is_variation_constrained_basis)

      if (is_basis_drop_possible) {
        # Remove basis functions whose coefficients need to be both positive
        # and negative
        is_nonzero <- !(is_positive_basis & is_negative_basis)
        basis_names <- basis_names[is_nonzero]
        basis_orders <- basis_orders[is_nonzero]
        is_positive_basis <- is_positive_basis[is_nonzero]
        is_negative_basis <- is_negative_basis[is_nonzero]
        is_variation_constrained_basis <- (
          is_variation_constrained_basis[is_nonzero]
        )
      }
    }

    basis_names <- append(basis_names, full_order_basis_names)
    is_positive_basis <- append(is_positive_basis, is_positive_full_order_basis)
    is_negative_basis <- append(is_negative_basis, is_negative_full_order_basis)
    is_variation_constrained_basis <- append(
      is_variation_constrained_basis, is_variation_constrained_full_order_basis
    )
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
  lasso_matrix <- lasso_matrix[, is_included_basis, drop = FALSE]
  is_positive_basis <- is_positive_basis[is_included_basis]
  is_negative_basis <- is_negative_basis[is_included_basis]
  is_variation_constrained_basis <- (
    is_variation_constrained_basis[is_included_basis]
  )

  # ============================================================================
  list(
    lasso_matrix = lasso_matrix,
    is_positive_basis = is_positive_basis,
    is_negative_basis = is_negative_basis,
    is_variation_constrained_basis = is_variation_constrained_basis,
    is_included_basis = is_included_basis
  )
}


get_lasso_matrix_emhk_nonlattice <- function(X_eval, X_design,
                                             max_vals, min_vals, s,
                                             increasing_covariates,
                                             decreasing_covariates,
                                             variation_constrained_covariates,
                                             is_included_basis = NULL) {
  d <- ncol(X_design)
  unique_entries <- get_unique_column_entries(X_design, 'emhk')
  for (col in (1L:d)) {
    if (length(unique_entries[[col]]) == 0L) {
      stop(paste0('All the values of "', colnames(X_design)[col], '" are zero. Please remove that variable.'))
    }
  }

  # Consider the constant term =================================================
  lasso_matrix <- matrix(rep(1.0, nrow(X_eval)), ncol = 1L)
  # Give a name to the basis function
  colnames(lasso_matrix) <- c("(Intercept)")
  # Record whether each basis function is restricted to have a positive
  # coefficient or a negative coefficient and whether its variation is
  # constrained or not
  is_positive_basis <- c(FALSE)
  is_negative_basis <- c(FALSE)
  is_variation_constrained_basis <- c(FALSE)
  # ============================================================================

  # Add the first order terms ==================================================
  for (col in (1L:d)) {
    column_unique <- unique_entries[[col]]

    X_eval_col <- X_eval[, col]
    basis <- sapply(column_unique, simplify = TRUE, function(entry) {
      compute_indicator(X_eval_col - entry)
    })
    if (is.vector(basis)) {
      basis <- matrix(basis, nrow = 1)
    }

    basis_names <- sapply(column_unique, simplify = TRUE, function(entry) {
      paste0("I(", colnames(X_design)[col], "-",
             scale_back_matrix_entry(entry, max_vals[col], min_vals[col],
                                     digits = 4L),
             ")")
    })
    colnames(basis) <- basis_names

    lasso_matrix <- cbind(lasso_matrix, basis)
    is_positive_basis <- c(
      is_positive_basis,
      rep((col %in% increasing_covariates), length(column_unique))
    )
    is_negative_basis <- c(
      is_negative_basis,
      rep((col %in% decreasing_covariates), length(column_unique))
    )
    is_variation_constrained_basis <- c(
      is_variation_constrained_basis,
      rep((col %in% variation_constrained_covariates), length(column_unique))
    )
  }
  # ============================================================================

  # Add the higher order terms =================================================
  if (s >= 2L) {
    for (order in (2L:s)) {
      subsets <- utils::combn(d, s)  # all subsets of {1, ... , d} of size s

      for (index in ncol(subsets)) {
        subset <- subsets[, index]  # a particular subset of size s

        # Check and record which constraints are induced by the selected
        # covariates
        is_positive <- any(subset %in% increasing_covariates)
        is_negative <- any(subset %in% decreasing_covariates)
        is_variation_constrained <- any(
          subset %in% variation_constrained_covariates
        )
        if (is_positive && is_negative) next

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
        X_eval_subset <- X_eval[, subset, drop = FALSE]
        basis <- sapply((1L:nrow(X_min)), simplify = TRUE, function(row) {
          X_min_row <- X_min[row, ]
          apply(X_eval_subset, MARGIN = 1L, FUN = function(eval_point) {
            all(compute_indicator(eval_point - X_min_row))
          })
        })
        if (is.vector(basis)) {
          basis <- matrix(basis, nrow = 1)
        }

        # Give names to the basis functions
        basis_names <- sapply((1L:nrow(X_min)), simplify = TRUE, function(row) {
          basis_name <- ""
          for (col in (1L:ncol(X_min))) {
            basis_name <- paste0(basis_name, "I(",
                                 colnames(X_design)[subset[col]], "-",
                                 scale_back_matrix_entry(
                                   X_min[row, col],
                                   max_vals[subset[col]], min_vals[subset[col]],
                                   digits = 4L
                                 ),
                                 ")")
          }
          basis_name
        })
        colnames(basis) <- basis_names

        lasso_matrix <- cbind(lasso_matrix, basis)
        is_positive_basis <- c(
          is_positive_basis, rep(is_positive, nrow(X_min))
        )
        is_negative_basis <- c(
          is_negative_basis, rep(is_negative, nrow(X_min))
        )
        is_variation_constrained_basis <- c(
          is_variation_constrained_basis,
          rep(is_variation_constrained, nrow(X_min))
        )
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
  lasso_matrix <- lasso_matrix[, is_included_basis, drop = FALSE]
  is_positive_basis <- is_positive_basis[is_included_basis]
  is_negative_basis <- is_negative_basis[is_included_basis]
  is_variation_constrained_basis <- (
    is_variation_constrained_basis[is_included_basis]
  )

  # ============================================================================
  list(
    lasso_matrix = lasso_matrix,
    is_positive_basis = is_positive_basis,
    is_negative_basis = is_negative_basis,
    is_variation_constrained_basis = is_variation_constrained_basis,
    is_included_basis = is_included_basis
  )
  # ============================================================================
}


get_lasso_matrix_tcmars <- function(X_eval, X_design, max_vals, min_vals, s,
                                    number_of_bins,
                                    concave_covariates, convex_covariates,
                                    variation_constrained_covariates,
                                    extra_linear_covariates,
                                    is_included_basis = NULL) {
  # Record whether basis dropping can happen due to the sign constraints
  is_basis_drop_possible <- !(is.null(concave_covariates)
                              || is.null(convex_covariates))

  # Extract the extra linear covariates from the matrix. They will be added back
  # at the end.
  if (!is.null(extra_linear_covariates)) {
    lasso_matrix_extra_linear <- X_eval[, extra_linear_covariates, drop = FALSE]
    colnames(lasso_matrix_extra_linear) <- (
      colnames(X_design)[extra_linear_covariates]
    )
    basis_scale_factors_extra_linear <- (
      max_vals[extra_linear_covariates] - min_vals[extra_linear_covariates]
    )

    X_eval <- X_eval[, -extra_linear_covariates, drop = FALSE]
    X_design <- X_design[, -extra_linear_covariates, drop = FALSE]
    max_vals <- max_vals[-extra_linear_covariates]
    min_vals <- min_vals[-extra_linear_covariates]
    number_of_bins <- number_of_bins[-extra_linear_covariates]
  }

  d <- ncol(X_design)
  unique_entries <- get_unique_column_entries(X_design, 'tcmars', number_of_bins)
  for (col in (1L:d)) {
    if (length(unique_entries[[col]]) == 0L) {
      stop(paste0('All the values of "', colnames(X_design)[col], '" are zero. Please remove that variable.'))
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

    if (is_basis_drop_possible) {
      # Record whether each hinge function is restricted to have a positive
      # coefficient or a negative coefficient and whether its variation is
      # constrained or not
      is_positive_hinge <- lapply((1L:d), function(col) {
        if (col %in% convex_covariates) {
          c(FALSE, FALSE, rep(TRUE, length(unique_entries[[col]]) - 1L))
        } else {
          rep(FALSE, length(unique_entries[[col]]) + 1L)
        }
      })
      is_negative_hinge <- lapply((1L:d), function(col) {
        if (col %in% concave_covariates) {
          c(FALSE, FALSE, rep(TRUE, length(unique_entries[[col]]) - 1L))
        } else {
          rep(FALSE, length(unique_entries[[col]]) + 1L)
        }
      })
      is_variation_constrained_hinge <- lapply((1L:d), function(col) {
        if (col %in% variation_constrained_covariates) {
          c(FALSE, FALSE, rep(TRUE, length(unique_entries[[col]]) - 1L))
        } else {
          rep(FALSE, length(unique_entries[[col]]) + 1L)
        }
      })
    }

    # Evaluate basis functions at the (row)-th evaluation point
    basis <- hinges[[1L]]
    # Compute the order of each basis function which is defined as the number of
    # hinge functions multiplied
    basis_orders <- hinges_orders[[1L]]
    if (is_basis_drop_possible) {
      # Record whether each basis function is restricted to have a positive
      # coefficient or a negative coefficient and whether its variation is
      # constrained or not
      is_positive_basis <- is_positive_hinge[[1L]]
      is_negative_basis <- is_negative_hinge[[1L]]
      is_variation_constrained_basis <- is_variation_constrained_hinge[[1L]]
    }

    # Put aside basis functions that reach the maximum allowed order s in the
    # for loop below to reduce computational cost
    full_order_basis <- c()

    if (d >= 2L) {
      for (k in (2L:d)) {
        # Put aside basis functions that reach the maximum allowed order s
        is_full_order <- (basis_orders == s)
        full_order_basis <- append(full_order_basis, basis[is_full_order])

        basis <- basis[!is_full_order]
        basis_orders <- basis_orders[!is_full_order]

        basis <- outer(basis, hinges[[k]], "*")
        basis <- c(basis)
        basis_orders <- outer(basis_orders, hinges_orders[[k]], "+")
        basis_orders <- c(basis_orders)

        if (is_basis_drop_possible) {
          is_positive_basis <- is_positive_basis[!is_full_order]
          is_negative_basis <- is_negative_basis[!is_full_order]
          is_variation_constrained_basis <- (
            is_variation_constrained_basis[!is_full_order]
          )

          # The descendants of basis functions that have positive (resp.
          # negative) coefficients are required to have positive (resp. negative)
          # coefficients. The descendants of basis functions whose variation is
          # constrained are also under the variation constraint
          is_positive_basis <- outer(is_positive_basis,
                                     is_positive_hinge[[k]], "|")
          is_positive_basis <- c(is_positive_basis)
          is_negative_basis <- outer(is_negative_basis,
                                     is_negative_hinge[[k]], "|")
          is_negative_basis <- c(is_negative_basis)
          is_variation_constrained_basis <- outer(
            is_variation_constrained_basis,
            is_variation_constrained_hinge[[k]], "|"
          )
          is_variation_constrained_basis <- c(is_variation_constrained_basis)

          # Remove basis functions whose coefficients need to be both positive
          # and negative
          is_nonzero <- !(is_positive_basis & is_negative_basis)
          basis <- basis[is_nonzero]
          basis_orders <- basis_orders[is_nonzero]
          is_positive_basis <- is_positive_basis[is_nonzero]
          is_negative_basis <- is_negative_basis[is_nonzero]
          is_variation_constrained_basis <- (
            is_variation_constrained_basis[is_nonzero]
          )
        }
      }
    }

    basis <- append(basis, full_order_basis)
  }))
  # ============================================================================

  # Create column names and find the constrained columns of the matrix =========
  # Give a name to each hinge function
  hinges_names <- lapply((1L:d), function(col) {
    names <- sapply((1L:length(unique_entries[[col]])), function(k) {
      paste0("H(", colnames(X_design)[col], "-",
             scale_back_matrix_entry(unique_entries[[col]][k],
                                     max_vals[col], min_vals[col],
                                     digits = 4L),
             ")")
      })
    c("", names)
  })

  # Record the order of each hinge function. We say constant functions have
  # order zero and hinge functions have order one.
  hinges_orders <- lapply((1L:d), function(col) {
    c(0L, rep(1L, length(unique_entries[[col]])))
  })

  # Record whether each hinge function is restricted to have a positive
  # coefficient or a negative coefficient and whether its variation is
  # constrained or not
  is_positive_hinge <- lapply((1L:d), function(col) {
    if (col %in% convex_covariates) {
      c(FALSE, FALSE, rep(TRUE, length(unique_entries[[col]]) - 1L))
    } else {
      rep(FALSE, length(unique_entries[[col]]) + 1L)
    }
  })
  is_negative_hinge <- lapply((1L:d), function(col) {
    if (col %in% concave_covariates) {
      c(FALSE, FALSE, rep(TRUE, length(unique_entries[[col]]) - 1L))
    } else {
      rep(FALSE, length(unique_entries[[col]]) + 1L)
    }
  })
  is_variation_constrained_hinge <- lapply((1L:d), function(col) {
    if (col %in% variation_constrained_covariates) {
      c(FALSE, FALSE, rep(TRUE, length(unique_entries[[col]]) - 1L))
    } else {
      rep(FALSE, length(unique_entries[[col]]) + 1L)
    }
  })

  # Record the scale factor for each hinge function
  hinges_scale_factors <- lapply((1L:d), function(col) {
    c(1.0, rep(1.0 / (max_vals[col] - min_vals[col]), length(unique_entries[[col]])))
  })

  # Give a name to each basis function
  basis_names <- hinges_names[[1L]]
  # Compute the order of each basis function which is defined as the number of
  # hinge functions multiplied
  basis_orders <- hinges_orders[[1L]]
  # Record whether each basis function is restricted to have a positive
  # coefficient or a negative coefficient and whether its variation is
  # constrained or not
  is_positive_basis <- is_positive_hinge[[1L]]
  is_negative_basis <- is_negative_hinge[[1L]]
  is_variation_constrained_basis <- is_variation_constrained_hinge[[1L]]
  # Compute the scale factor for each basis function
  basis_scale_factors <- hinges_scale_factors[[1L]]

  # Put aside basis functions that reach the maximum allowed order s in the
  # for loop below to reduce computational cost
  full_order_basis_names <- c()
  is_positive_full_order_basis <- c()
  is_negative_full_order_basis <- c()
  is_variation_constrained_full_order_basis <- c()
  full_order_basis_scale_factors <- c()

  if (d >= 2L) {
    for (k in (2L:d)) {
      # Put aside basis functions that reach the maximum allowed order s
      is_full_order <- (basis_orders == s)
      full_order_basis_names <- append(
        full_order_basis_names, basis_names[is_full_order]
      )
      is_positive_full_order_basis <- append(
        is_positive_full_order_basis, is_positive_basis[is_full_order]
      )
      is_negative_full_order_basis <- append(
        is_negative_full_order_basis, is_negative_basis[is_full_order]
      )
      is_variation_constrained_full_order_basis <- append(
        is_variation_constrained_full_order_basis,
        is_variation_constrained_basis[is_full_order]
      )
      full_order_basis_scale_factors <- append(
        full_order_basis_scale_factors, basis_scale_factors[is_full_order]
      )

      basis_names <- basis_names[!is_full_order]
      basis_orders <- basis_orders[!is_full_order]
      is_positive_basis <- is_positive_basis[!is_full_order]
      is_negative_basis <- is_negative_basis[!is_full_order]
      is_variation_constrained_basis <- (
        is_variation_constrained_basis[!is_full_order]
      )
      basis_scale_factors <- basis_scale_factors[!is_full_order]

      basis_names <- outer(basis_names, hinges_names[[k]], "paste0")
      basis_names <- c(basis_names)
      basis_orders <- outer(basis_orders, hinges_orders[[k]], "+")
      basis_orders <- c(basis_orders)

      # The descendants of basis functions that have positive (resp.
      # negative) coefficients are required to have positive (resp. negative)
      # coefficients. The descendants of basis functions whose variation is
      # constrained are also under the variation constraint
      is_positive_basis <- outer(is_positive_basis, is_positive_hinge[[k]], "|")
      is_positive_basis <- c(is_positive_basis)
      is_negative_basis <- outer(is_negative_basis, is_negative_hinge[[k]], "|")
      is_negative_basis <- c(is_negative_basis)
      is_variation_constrained_basis <- outer(
        is_variation_constrained_basis, is_variation_constrained_hinge[[k]], "|"
      )
      is_variation_constrained_basis <- c(is_variation_constrained_basis)
      basis_scale_factors <- outer(
        basis_scale_factors, hinges_scale_factors[[k]], "*"
      )
      basis_scale_factors <- c(basis_scale_factors)

      if (is_basis_drop_possible) {
        # Remove basis functions whose coefficients need to be both positive
        # and negative
        is_nonzero <- !(is_negative_basis & is_positive_basis)
        basis_names <- basis_names[is_nonzero]
        basis_orders <- basis_orders[is_nonzero]
        is_positive_basis <- is_positive_basis[is_nonzero]
        is_negative_basis <- is_negative_basis[is_nonzero]
        is_variation_constrained_basis <- (
          is_variation_constrained_basis[is_nonzero]
        )
        basis_scale_factors <- basis_scale_factors[is_nonzero]
      }
    }

    basis_names <- append(basis_names, full_order_basis_names)
    is_positive_basis <- append(is_positive_basis, is_positive_full_order_basis)
    is_negative_basis <- append(is_negative_basis, is_negative_full_order_basis)
    is_variation_constrained_basis <- append(
      is_variation_constrained_basis, is_variation_constrained_full_order_basis
    )
    basis_scale_factors <- append(
      basis_scale_factors, full_order_basis_scale_factors
    )
  }

  basis_names[1L] <- "(Intercept)"
  colnames(lasso_matrix) <- basis_names

  # Add back the extra linear covariates to the matrix
  if (!is.null(extra_linear_covariates)) {
    lasso_matrix <- as.matrix(cbind(lasso_matrix, lasso_matrix_extra_linear))
    is_positive_basis <- c(
      is_positive_basis, rep(FALSE, length(extra_linear_covariates))
    )
    is_negative_basis <- c(
      is_negative_basis, rep(FALSE, length(extra_linear_covariates))
    )
    is_variation_constrained_basis <- c(
      is_variation_constrained_basis, rep(FALSE, length(extra_linear_covariates))
    )
    basis_scale_factors <- c(
      basis_scale_factors, basis_scale_factors_extra_linear
    )
  }

  if (is.null(is_included_basis)) {
    # Remove all zero columns
    is_included_basis <- apply((lasso_matrix != 0), MARGIN = 2, any)
  } else {
    if (length(is_included_basis) != ncol(lasso_matrix)) {
      stop('`length(is_included_basis)` must be equal to the number of columns of the LASSO matrix.')
    }
  }

  # Include specified columns
  lasso_matrix <- lasso_matrix[, is_included_basis, drop = FALSE]
  is_positive_basis <- is_positive_basis[is_included_basis]
  is_negative_basis <- is_negative_basis[is_included_basis]
  is_variation_constrained_basis <- (
    is_variation_constrained_basis[is_included_basis]
  )
  basis_scale_factors <- basis_scale_factors[is_included_basis]

  # ============================================================================
  list(
    lasso_matrix = lasso_matrix,
    is_positive_basis = is_positive_basis,
    is_negative_basis = is_negative_basis,
    is_variation_constrained_basis = is_variation_constrained_basis,
    basis_scale_factors = basis_scale_factors,
    is_included_basis = is_included_basis
  )
}


scale_back_matrix_entry <- function(entry, max_val, min_val, digits = 4L) {
  format((max_val - min_val) * entry + min_val, digits = digits)
}
