#' Fit to data a nonparametric regression model with mixed derivative constraints
#'
#' Given an estimation method, this function builds a model by solving the
#' corresponding constrained LASSO problem. Available estimation methods are
#' entirely monotonic regression ("em"), Hardy—Krause variation denoising ("hk"),
#' their generalization ("emhk"), totally concave regression ("tc"), MARS via
#' LASSO ("mars"), and their generalization ("tcmars"). For details about the
#' corresponding LASSO problems, see, for example, Section 3 of Fang et al.
#' (2021) (for entirely monotonic regression and Hardy—Krause variation
#' denoising) and Section 2 of Ki et al. (2024+) (for MARS via LASSO). There are
#' some other ongoing research and working papers, and they will be available in
#' the future.
#'
#' @param X_design A numeric design matrix. Each row corresponds to an
#'   individual data.
#' @param y A numeric observation vector of a response variable.
#' @param s A numeric scalar indicating the maximum order of interaction between
#'   covariates allowed in the estimation method.
#' @param method A string indicating the estimation method. One of "em", "hk",
#'   "emhk", "tc", "mars", and "tcmars".
#' @param V A numeric scalar. An upper bound on complexity measure (variation)
#'   in a scaled domain. Required for "hk" and "mars", and possibly used for
#'   "emhk" and "tcmars".
#' @param threshold A numeric scalar to determine whether each component of the
#'   solution to the LASSO problem is zero or not.
#' @param is_scaled A logical scalar for whether the design matrix is scaled so
#'   that every entry is between 0 and 1. If `FALSE`, the min-max scaling is
#'   applied to each column of the design matrix.
#' @param is_lattice A logical scalar for whether the design is a lattice or not.
#'   Only used for "em", "hk", and "emhk". See details below.
#' @param is_monotonically_increasing A logical scalar for whether the method is
#'   entirely monotonically increasing regression or entirely monotonically
#'   decreasing regression. Only used for "em".
#' @param is_totally_concave A logical scalar for whether the method is totally
#'   concave regression or totally convex regression. Only used for "tc".
#' @param number_of_bins An integer or an integer vector of the numbers of bins
#'   for the approximate method. Currently available for "tc", "mars", and
#'   "tcmars". See details below.
#' @param increasing_covariates An integer or a string vector of monotonically
#'   increasing covariates. Possibly used for "em" and "emhk". See details below.
#' @param decreasing_covariates An integer or a string vector of monotonically
#'   decreasing covariates. Possibly used for "em" and "emhk". See details below.
#' @param concave_covariates An integer or a string vector of concave covariates.
#'   Possibly used for "tc" and "tcmars". See details below.
#' @param convex_covariates An integer or a string vector of convex covariates.
#'   Possibly used for "tc" and "tcmars". See details below.
#' @param extra_linear_covariates An integer vector or a string vector of extra
#'   linear covariates added to the model. Possibly used for "tc", "mars", and
#'   "tcmars". See details below.
#' @details
#' You can also set `is_lattice` as `TRUE` even when the design is not a lattice.
#' If `is_lattice` is `TRUE`, the model is fitted in a faster way; whereas if it
#' is `FALSE`, the model is fitted in a memory-efficient way. Currently,
#' `is_lattice` can be used for entirely monotonic regression, Hardy—Krause
#' variation denoising, and their generalization.
#'
#' The approximate method is used if `number_of_bins` is not `NULL`. Currently,
#' the approximate method is available for totally concave regression, MARS via
#' LASSO, and their generalization. You can put an integer into `number_of_bins`
#' if the numbers of bins for the approximate method are the same for all
#' covariates. If they are different, you need to put an integer vector of the
#' numbers of bins. The number of bins can be `NA` for some covariates. In that
#' case, the approximate method is not applied to such covariates.
#'
#' Using `extra_linear_covariates`, you can add to the model extra covariates
#' that have linear relationship with the response variable. The interaction
#' between these covariates and the interaction between these covariates and the
#' other covariates are not considered in the model.
#'
#' For entirely monotonic regression (resp. totally concave regression), you can
#' fit a model using both covariates that have a increasing (resp. concave)
#' relationship with the response variable and covariates that have a decreasing
#' (resp. convex) relationship with the response variable. You can utilize
#' `increasing_covariates` and `decreasing_covariates` for entirely monotonic
#' regression, and `concave_covariates` and `convex_covariates` for totally
#' concave regression.
#'
#' For entirely monotonic regression, if exactly one of `increasing_covariates`
#' and `decreasing_covariates` is not `NULL`, the one that is `NULL` is
#' automatically changed to the complement of the other. If both of them are not
#' `NULL`, their union must be the whole set of covariates.
#'
#' For totally concave regression, if exactly one of `concave_covariates` and
#' `convex_covariates` is not `NULL`, the one that is `NULL` is automatically
#' converted to the relative complement of the other in the complement of
#' `extra_linear_covariates`. If both of them are not `NULL`, their union must
#' include all covariates except those in `extra_linear_covariates`.
#'
#' In the generalization of entirely monotonic regression and Hardy—Krause
#' variation denoising, you can also have covariates whose variation is
#' constrained as in Hardy—Krause variation denoising. The covariates that are
#' neither in `increasing_covariates` nor `decreasing_covariates`, and all
#' interactions including at least one of them are under the variation constraint.
#'
#' Similarly, in the generalization of totally concave regression and MARS via
#' LASSO, you can also have covariates whose variation is constrained as in
#' MARS via LASSO. The variation constraint is imposed to the covariates that do
#' not belong to any of `concave_covariates`, `convex_covariates`, and
#' `extra_linear_covariates`, and all interactions including at least one of them.
#'
#' For `increasing_covariates`, `decreasing_covariates`, `concave_covariates`,
#' `convex_covariates`, and `extra_linear_covariates`, you can put either an
#' integer vector of column indices or a string vector of column names.
#'
#' @references Ki, D., Fang, B., and Guntuboyina, A. (2024+). MARS via LASSO.
#'   Accepted at \emph{Annals of Statistics}. Available at
#'   \url{https://arxiv.org/abs/2111.11694}.
#' @references Fang, B., Guntuboyina, A., and Sen, B. (2021). Multivariate
#'   extensions of isotonic regression and total variation denoising via entire
#'   monotonicity and Hardy—Krause variation. \emph{Annals of Statistics},
#'   \strong{49}(2), 769-792.
#' @examples
#' fstar <- function(x) {(
#'   (x[1] - 0.25 >= 0) + (x[2] - 0.25 >= 0)
#'   + (x[1] - 0.25 >= 0) * (x[2] - 0.25 >= 0)
#' )}  # the true underlying function
#' X_design <- expand.grid(rep(list(seq(0, 1, length.out = 5L)), 3L))
#' colnames(X_design) <- c("VarA", "VarB", "VarC")
#' theta <- apply(X_design, MARGIN = 1L, FUN = fstar)
#' sigma <- 0.1
#' y <- theta + sigma * rnorm(nrow(X_design))
#'
#' regmdc(X_design, y, s = 1L, method = "em")
#' regmdc(X_design, y, s = 2L, method = "em")
#' regmdc(X_design, y, s = 2L, method = "em", is_monotonically_increasing = FALSE)
#' regmdc(X_design, y, s = 2L, method = "em", increasing_covariates = c(1L, 2L))
#' regmdc(X_design, y, s = 2L, method = "em", decreasing_covariates = c(3L))
#' regmdc(X_design, y, s = 2L, method = "em",
#'        increasing_covariates = c("VarA", "VarB"),
#'        decreasing_covariates = c("VarC"))
#' regmdc(X_design, y, s = 2L, method = "em", is_scaled = TRUE)
#' regmdc(X_design, y, s = 2L, method = "em", is_lattice = TRUE)
#' regmdc(X_design, y, s = 2L, method = "hk", V = 3.0)
#' regmdc(X_design, y, s = 2L, method = "emhk", V = 2.0,
#'        increasing_covariates = c(1L))
#' regmdc(X_design, y, s = 2L, method = "emhk", V = 3.0,
#'        decreasing_covariates = c(3L))
#' regmdc(X_design, y, s = 2L, method = "emhk", V = 2.0,
#'        increasing_covariates = c("VarA"), decreasing_covariates = c("VarC"))
#'
#' fstar <- function(x) {(
#'   - max(x[1] - 0.25, 0) - max(x[2] - 0.25, 0)
#'   - max(x[1] - 0.25, 0) * max(x[2] - 0.25, 0)
#' )}  # the true underlying function
#' X_design <- cbind(runif(50), runif(50), runif(50))
#' colnames(X_design) <- c("VarA", "VarB", "VarC")
#' theta <- apply(X_design, MARGIN = 1L, FUN = fstar)
#' sigma <- 0.1
#' y <- theta + sigma * rnorm(nrow(X_design))
#'
#' regmdc(X_design, y, s = 2L, method = "tc")
#' regmdc(X_design, y, s = 2L, method = "tc", is_totally_concave = FALSE)
#' regmdc(X_design, y, s = 2L, method = "tc", concave_covariates = c(1L, 2L))
#' regmdc(X_design, y, s = 2L, method = "tc", convex_covariates = c(3L))
#' regmdc(X_design, y, s = 2L, method = "tc",
#'        concave_covariates = c("VarA", "VarB"), convex_covariates = c("VarC"))
#' regmdc(X_design, y, s = 2L, method = "tc", extra_linear_covariates = c(3L))
#' regmdc(X_design, y, s = 2L, method = "tc", is_totally_concave = FALSE,
#'        extra_linear_covariates = c("VarC"))
#' regmdc(X_design, y, s = 2L, method = "tc", extra_linear_covariates = c(2L, 3L))
#' regmdc(X_design, y, s = 2L, method = "tc", concave_covariates = c("VarA"),
#'        extra_linear_covariates = c("VarC"))
#' regmdc(X_design, y, s = 2L, method = "tc", number_of_bins = 20L,
#'        extra_linear_covariates = c(3L))
#' regmdc(X_design, y, s = 2L, method = "mars", V = 3.0)
#' regmdc(X_design, y, s = 2L, method = "mars", V = 3.0, number_of_bins = 20L)
#' regmdc(X_design, y, s = 2L, method = "mars", V = 3.0,
#'        number_of_bins = c(10L, 20L, 20L))
#' regmdc(X_design, y, s = 2L, method = "mars", V = 3.0,
#'        number_of_bins = c(10L, 20L, NA))
#' regmdc(X_design, y, s = 2L, method = "mars", V = 3.0,
#'        number_of_bins = c(10L, 20L, NA), extra_linear_covariates = c("VarC"))
#' regmdc(X_design, y, s = 2L, method = "tcmars", V = 2.0,
#'        concave_covariates = c(1L))
#' regmdc(X_design, y, s = 2L, method = "tcmars", V = 2.0,
#'        concave_covariates = c(1L), convex_covariates = c(3L))
#' regmdc(X_design, y, s = 2L, method = "tcmars", V = 2.0,
#'        concave_covariates = c("VarA"), extra_linear_covariates = c("VarC"))
#' regmdc(X_design, y, s = 2L, method = "tcmars", number_of_bins = 20L,
#'        concave_covariates = c("VarA"), extra_linear_covariates = c("VarC"))
#' @export
regmdc <- function(X_design, y, s, method, V = Inf, threshold = 1e-6,
                   is_scaled = FALSE, is_lattice = FALSE,
                   is_monotonically_increasing = TRUE, is_totally_concave = TRUE,
                   number_of_bins = NULL,
                   increasing_covariates = NULL, decreasing_covariates = NULL,
                   concave_covariates = NULL, convex_covariates = NULL,
                   extra_linear_covariates = NULL) {
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
  if (!is.null(colnames(X_design))) {
    if (any(duplicated(colnames(X_design)))) {
      stop('The column names of `X_design` must be different.')
    }
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

  if (!is.logical(is_scaled)) {
    stop('`is_scaled` must be TRUE or FALSE.')
  }
  if (is_scaled) {
    if (min(X_design) < 0 || max(X_design) > 1) {
      stop('If `is_scaled` is TRUE, all the entries of `X_design` must be between 0 and 1.')
    }
  }

  if (!is.logical(is_lattice)) {
    stop('`is_lattice` must be TRUE or FALSE.')
  }

  if (!is.logical(is_monotonically_increasing)) {
    stop('`is_monotonically_increasing` must be TRUE or FALSE.')
  }

  if (!is.logical(is_totally_concave)) {
    stop('`is_totally_concave` must be TRUE or FALSE.')
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

  # If `increasing_covariates` is a string vector of column names, covert it
  # into an integer vector of column indices
  if (!is.null(increasing_covariates)) {
    if (method %in% c('tc', 'mars', 'tcmars', 'hk')) {
      stop('`increasing_covariates` can only be used for `em` and `emhk`.')
    }

    increasing_covariates <- unique(increasing_covariates)
    if (is.numeric(increasing_covariates)) {
      if (!is.integer(increasing_covariates)) {
        stop('`increasing_covariates` must be an integer vector.')
      }

      increasing_covariates <- sort(increasing_covariates)
      if ((increasing_covariates[1] <= 0)
          | (increasing_covariates[length(increasing_covariates)] > ncol(X_design))) {
        stop('Each integer in `increasing_covariates` must be at least 1 and at most `ncol(X_design)`.')
      }
    } else {
      increasing_covariates <- sapply(increasing_covariates, function(col_name) {
        col_name_index <- which(colnames(X_design) == col_name)
        if (length(col_name_index) == 0) {
          stop(paste0('`X_design` does not have a column with the name "', col_name, '".'))
        } else {
          col_name_index
        }
      })
    }
  }

  # If `decreasing_covariates` is a string vector of column names, covert it
  # into an integer vector of column indices
  if (!is.null(decreasing_covariates)) {
    if (method %in% c('tc', 'mars', 'tcmars', 'hk')) {
      stop('`decreasing_covariates` can only be used for `em` and `emhk`.')
    }

    decreasing_covariates <- unique(decreasing_covariates)
    if (is.numeric(decreasing_covariates)) {
      if (!is.integer(decreasing_covariates)) {
        stop('`decreasing_covariates` must be an integer vector.')
      }

      decreasing_covariates <- sort(decreasing_covariates)
      if ((decreasing_covariates[1] <= 0)
          | (decreasing_covariates[length(decreasing_covariates)] > ncol(X_design))) {
        stop('Each integer in `decreasing_covariates` must be at least 1 and at most `ncol(X_design)`.')
      }
    } else {
      decreasing_covariates <- sapply(decreasing_covariates, function(col_name) {
        col_name_index <- which(colnames(X_design) == col_name)
        if (length(col_name_index) == 0) {
          stop(paste0('`X_design` does not have a column with the name "', col_name, '".'))
        } else {
          col_name_index
        }
      })
    }
  }

  # Specify `increasing_covariates`, `decreasing_covariates`, and
  # `variation_constrained_covariates` when `method` is "em", "hk", or "emhk"
  if (method %in% c("em", "hk", "emhk")) {
    if (method == "em") {
      if (is.null(increasing_covariates) && is.null(decreasing_covariates)) {
        if (is_monotonically_increasing) {
          increasing_covariates <- (1L:ncol(X_design))
          decreasing_covariates <- NULL
        } else {
          increasing_covariates <- NULL
          decreasing_covariates <- (1L:ncol(X_design))
        }
      } else if (is.null(decreasing_covariates)) {
        decreasing_covariates <- (1L:ncol(X_design))[-increasing_covariates]
        if (length(increasing_covariates) == 0) {
          increasing_covariates <- NULL
        }
      } else if (is.null(increasing_covariates)){
        increasing_covariates <- (1L:ncol(X_design))[-decreasing_covariates]
        if (length(decreasing_covariates) == 0) {
          decreasing_covariates <- NULL
        }
      } else {
        if (length(unique(c(increasing_covariates, decreasing_covariates))) < ncol(X_design)) {
          stop('If `method` = `em` and neither `increasing_covariates` nor `decreasing_covariates` is NULL, all covariates must be included in at least one of them.')
        }
      }
      variation_constrained_covariates <- NULL
    } else if (method == "hk") {
      increasing_covariates <- NULL
      decreasing_covariates <- NULL
      variation_constrained_covariates <- (1L:ncol(X_design))
    } else {
      if (is.null(increasing_covariates) && is.null(decreasing_covariates)) {
        variation_constrained_covariates <- (1L:ncol(X_design))
      } else {
        variation_constrained_covariates <- (
          (1L:ncol(X_design))[-c(increasing_covariates, decreasing_covariates)]
        )
      }
    }
  }

  # If `concave_covariates` is a string vector of column names, covert it into
  # an integer vector of column indices
  if (!is.null(concave_covariates)) {
    if (method %in% c('em', 'hk', 'emhk', 'mars')) {
      stop('`concave_covariates` can only be used for `tc` and `tcmars`.')
    }

    concave_covariates <- unique(concave_covariates)
    if (is.numeric(concave_covariates)) {
      if (!is.integer(concave_covariates)) {
        stop('`concave_covariates` must be an integer vector.')
      }

      concave_covariates <- sort(concave_covariates)
      if ((concave_covariates[1] <= 0)
          | (concave_covariates[length(concave_covariates)] > ncol(X_design))) {
        stop('Each integer in `concave_covariates` must be at least 1 and at most `ncol(X_design)`.')
      }
    } else {
      concave_covariates <- sapply(concave_covariates, function(col_name) {
        col_name_index <- which(colnames(X_design) == col_name)
        if (length(col_name_index) == 0) {
          stop(paste0('`X_design` does not have a column with the name "', col_name, '".'))
        } else {
          col_name_index
        }
      })
    }
  }

  # If `convex_covariates` is a string vector of column names, covert it into an
  # integer vector of column indices
  if (!is.null(convex_covariates)) {
    if (method %in% c('em', 'hk', 'emhk', 'mars')) {
      stop('`convex_covariates` can only be used for `tc` and `tcmars`.')
    }

    convex_covariates <- unique(convex_covariates)
    if (is.numeric(convex_covariates)) {
      if (!is.integer(convex_covariates)) {
        stop('`convex_covariates` must be an integer vector.')
      }

      convex_covariates <- sort(convex_covariates)
      if ((convex_covariates[1] <= 0)
          | (convex_covariates[length(convex_covariates)] > ncol(X_design))) {
        stop('Each integer in `convex_covariates` must be at least 1 and at most `ncol(X_design)`.')
      }
    } else {
      convex_covariates <- sapply(convex_covariates, function(col_name) {
        col_name_index <- which(colnames(X_design) == col_name)
        if (length(col_name_index) == 0) {
          stop(paste0('`X_design` does not have a column with the name "', col_name, '".'))
        } else {
          col_name_index
        }
      })
    }
  }

  # If `extra_linear_covariates` is a string vector of column names, covert it
  # into an integer vector of column indices
  if (!is.null(extra_linear_covariates)) {
    if (method %in% c('em', 'hk', 'emhk')) {
      stop('`extra_linear_covariates` can only be used for `tc`, `mars`, and `tcmars`.')
    }

    extra_linear_covariates <- unique(extra_linear_covariates)
    if (is.numeric(extra_linear_covariates)) {
      if (!is.integer(extra_linear_covariates)) {
        stop('`extra_linear_covariates` must be an integer vector.')
      }

      extra_linear_covariates <- sort(extra_linear_covariates)
      if ((extra_linear_covariates[1] <= 0)
          | (extra_linear_covariates[length(extra_linear_covariates)] > ncol(X_design))) {
        stop('Each integer in `extra_linear_covariates` must be at least 1 and at most `ncol(X_design)`.')
      }
    } else {
      extra_linear_covariates <- sapply(extra_linear_covariates, function(col_name) {
        col_name_index <- which(colnames(X_design) == col_name)
        if (length(col_name_index) == 0) {
          stop(paste0('`X_design` does not have a column with the name "', col_name, '".'))
        } else {
          col_name_index
        }
      })
    }

    if (length(extra_linear_covariates) == ncol(X_design)) {
      stop('There must be at least one covariate that is not in `extra_linear_covariates`.')
    }
  }

  if (any(extra_linear_covariates %in% concave_covariates)) {
    stop('`concave_covariates` and `extra_linear_covariates` must not have intersection.')
  }
  if (any(extra_linear_covariates %in% convex_covariates)) {
    stop('`convex_covariates` and `extra_linear_covariates` must not have intersection.')
  }

  # Specify `concave_covariates`, `convex_covariates`, and
  # `variation_constrained_covariates` when `method` is "tc", "mars", or "tcmars"
  if (method %in% c("tc", "mars", "tcmars")) {
    if (method == "tc") {
      if (is.null(concave_covariates) && is.null(convex_covariates)) {
        if (is_totally_concave) {
          if (is.null(extra_linear_covariates)) {
            concave_covariates <- (1L:ncol(X_design))
          } else {
            concave_covariates <- (1L:ncol(X_design))[-extra_linear_covariates]
          }
          convex_covariates <- NULL
        } else {
          concave_covariates <- NULL
          if (is.null(extra_linear_covariates)) {
            convex_covariates <- (1L:ncol(X_design))
          } else {
            convex_covariates <- (1L:ncol(X_design))[-extra_linear_covariates]
          }
        }
      } else if (is.null(convex_covariates)) {
        convex_covariates <- (
          (1L:ncol(X_design))[-c(concave_covariates, extra_linear_covariates)]
        )
        if (length(convex_covariates) == 0) {
          convex_covariates <- NULL
        }
      } else if (is.null(concave_covariates)){
        concave_covariates <- (
          (1L:ncol(X_design))[-c(convex_covariates, extra_linear_covariates)]
        )
        if (length(concave_covariates) == 0) {
          concave_covariates <- NULL
        }
      } else {
        if (length(unique(c(concave_covariates, convex_covariates)))
            < (ncol(X_design) - length(extra_linear_covariates))) {
          stop('If `method` = `tc` and neither `concave_covariates` nor `convex_covariates` is NULL, all covariates except those in `extra_linear_covariates` must be included in at least one of them.')
        }
      }
      variation_constrained_covariates <- NULL
    } else if (method == "mars") {
      concave_covariates <- NULL
      convex_covariates <- NULL
      if (is.null(extra_linear_covariates)) {
        variation_constrained_covariates <- (1L:ncol(X_design))
      } else {
        variation_constrained_covariates <- (
          (1L:ncol(X_design))[-extra_linear_covariates]
        )
      }
    } else {
      if (is.null(extra_linear_covariates) && is.null(concave_covariates) && is.null(convex_covariates)) {
        variation_constrained_covariates <- (1L:ncol(X_design))
      } else {
        variation_constrained_covariates <- (
          (1L:ncol(X_design))[-c(concave_covariates,
                                 convex_covariates,
                                 extra_linear_covariates)]
        )
      }
    }
  }

  # ============================================================================
  # Give names to the columns of the design matrix if there aren't
  if (is.null(colnames(X_design))) {
    colnames(X_design) <- paste0("Var", (1L:ncol(X_design)))
  }

  # Find the maximal and minimal values of each covariate
  if (is_scaled) {
    max_vals <- rep(1, ncol(X_design))
    min_vals <- rep(0, ncol(X_design))
  } else {
    max_vals <- apply(X_design, MARGIN = 2L, max)
    min_vals <- apply(X_design, MARGIN = 2L, min)
  }

  for (col in (1L:ncol(X_design))) {
    if (max_vals[col] == min_vals[col]) {
      stop(paste0('All the values of "', colnames(X_design)[col], '" are the same. Please remove that variable.'))
    }
  }

  # Obtain the matrix for the LASSO problem and the logical vectors for whether
  # each basis function is restricted to have a positive coefficient or a
  # negative coefficient and whether its variation is constrained or not in the
  # model. For totally concave regression, MARS via LASSO, and their
  # generalization, the scale factors of basis functions, which are needed for
  # rescaling, are additionally collected.
  matrix_with_additional_info <- get_lasso_matrix(
    X_design, X_design, max_vals, min_vals, s, method, is_scaled, is_lattice,
    number_of_bins, increasing_covariates, decreasing_covariates,
    concave_covariates, convex_covariates, variation_constrained_covariates,
    extra_linear_covariates
  )
  M <- matrix_with_additional_info$lasso_matrix
  is_positive_basis <- matrix_with_additional_info$is_positive_basis
  is_negative_basis <- matrix_with_additional_info$is_negative_basis
  is_variation_constrained_basis <- (
    matrix_with_additional_info$is_variation_constrained_basis
  )
  if (method %in% c('tc', 'mars', 'tcmars')) {
    basis_scale_factors <- matrix_with_additional_info$basis_scale_factors
  }
  is_included_basis <- matrix_with_additional_info$is_included_basis

  # Solve the LASSO problem
  solution <- solve_constrained_lasso(
    y, M, V = V,
    is_sum_constrained_component = is_variation_constrained_basis,
    is_positive_component = is_positive_basis,
    is_negative_component = is_negative_basis
  )

  names(solution) <- colnames(M)

  # Remove the zero components from the solution
  is_nonzero_component <- (abs(solution) >= threshold)
  coefficients <- solution[is_nonzero_component]

  # Find the components whose sum of absolute values is constrained, whose signs
  # are constrained, and which are unconstrained
  is_variation_constrained_basis <- (
    is_variation_constrained_basis[is_nonzero_component]
  )
  is_sign_constrained_basis <- (is_positive_basis | is_negative_basis)
  is_sign_constrained_basis <- (
    is_sign_constrained_basis[is_nonzero_component]
  )
  is_constrained_basis <- (is_variation_constrained_basis
                           | is_sign_constrained_basis)

  # Compute the variation of the fitted function in a scaled domain
  if (any(is_variation_constrained_basis)) {
    V_solution <- sum(abs(coefficients[is_variation_constrained_basis]))
  } else {
    V_solution <- NULL
  }

  # Rescale the coefficients if necessary
  if (method %in% c('em', 'hk', 'emhk')) {
    coefficients_rescaled <- coefficients
  } else {
    basis_scale_factors <- basis_scale_factors[is_nonzero_component]
    coefficients_rescaled <- coefficients * basis_scale_factors
  }

  variation_constrained_components <- (
    coefficients_rescaled[is_variation_constrained_basis]
  )
  if (length(variation_constrained_components) == 0) {
    variation_constrained_components <- NULL
  }
  sign_constrained_components <- coefficients_rescaled[is_sign_constrained_basis]
  if (length(sign_constrained_components) == 0) {
    sign_constrained_components <- NULL
  }
  unconstrained_components <- coefficients_rescaled[!is_constrained_basis]
  if (length(unconstrained_components) == 0) {
    unconstrained_components <- NULL
  }

  # Compute the fitted values at the design points
  fitted_values <- M[, is_nonzero_component, drop = FALSE] %*% coefficients

  # ============================================================================
  if (!is.null(increasing_covariates)) {
    names(increasing_covariates) <- colnames(X_design)[increasing_covariates]
  }
  if (!is.null(decreasing_covariates)) {
    names(decreasing_covariates) <- colnames(X_design)[decreasing_covariates]
  }
  if (!is.null(concave_covariates)) {
    names(concave_covariates) <- colnames(X_design)[concave_covariates]
  }
  if (!is.null(convex_covariates)) {
    names(convex_covariates) <- colnames(X_design)[convex_covariates]
  }
  if (!is.null(variation_constrained_covariates)) {
    names(variation_constrained_covariates) <- (
      colnames(X_design)[variation_constrained_covariates]
    )
  }
  if (!is.null(extra_linear_covariates)) {
    names(extra_linear_covariates) <- colnames(X_design)[extra_linear_covariates]
  }

  regmdc_model <- list(
    X_design = X_design,
    y = y,
    s = s,
    method = method,
    V = V,
    threshold = threshold,
    is_scaled = is_scaled,
    is_lattice = is_lattice,
    number_of_bins = number_of_bins,
    increasing_covariates = increasing_covariates,
    decreasing_covariates = decreasing_covariates,
    concave_covariates = concave_covariates,
    convex_covariates = convex_covariates,
    variation_constrained_covariates = variation_constrained_covariates,
    extra_linear_covariates = extra_linear_covariates,
    max_vals = max_vals,
    min_vals = min_vals,
    coefficients = coefficients,
    coefficients_rescaled = coefficients_rescaled,
    variation_constrained_components = variation_constrained_components,
    sign_constrained_components = sign_constrained_components,
    unconstrained_components = unconstrained_components,
    V_solution = V_solution,
    fitted_values = fitted_values,
    is_included_basis = is_included_basis,
    is_nonzero_component = is_nonzero_component
  )
  class(regmdc_model) <- "regmdc"

  regmdc_model
}
