#' Find the solution to the constrained LASSO problem specified below
#'
#' Using the conic optimization tools of \code{\link[Rmosek]{mosek}}, this
#' function solves the problem: min_x \eqn{|| y - M x ||^2} subject to
#' * \eqn{\sum_i |x_i| \le V} where the sum is over the components of x where
#'   `is_sum_constrained_component` is `TRUE`,
#' * the components of x where `is_positive_component` is `TRUE` are restricted
#'   to be positive
#' * the components of x where `is_negative_component` is `TRUE` are restricted
#'   to be negative.
#'
#' @param y A numeric vector.
#' @param M A numeric matrix.
#' @param V A numeric scalar.
#' @param is_sum_constrained_component A logical vector for whether each
#'   component of x appears in the constraint on the sum of absolute values.
#' @param is_positive_component A logical vector for whether each component of x
#'   is restricted to be positive.
#' @param is_negative_component A logical vector for whether each component of x
#'   is restricted to be negative.
#' @seealso \url{https://docs.mosek.com/9.3/rmosek/tutorial-cqo-shared.html}
#'   for more details about Mosek's conic optimization.
solve_constrained_lasso <- function(y, M, V = Inf,
                                    is_sum_constrained_component = NULL,
                                    is_positive_component = NULL,
                                    is_negative_component = NULL) {
  n <- dim(M)[1]
  p <- dim(M)[2]
  prob <- list(sense = "min")

  # Represent our optimization variable as the concatenation of the following:
  # * t (the scalar variable in the cone constraint),
  # * u (a vector of length n), to be set equal to y - Mx for our use in the
  #   cone constraint,
  # * x_+ (the positive parts of x) (a vector of length p),
  # * x_- (the negative parts of x) (a vector of length p).

  # Set the scalar t in the cone constraint as our objective function
  prob$c <- c(1, rep(0, n + 2L * p))

  # Define the cone constraint $t \ge \sum_i u_i^2$
  prob$cones <- matrix(list(), nrow = 2L, ncol = 1L)
  rownames(prob$cones) <- c("type", "sub")
  prob$cones[, 1] <- list("QUAD", c(1L, 2L:(n + 1L)))

  # Define the sign constraints of x. The positive and the negative parts of x
  # need to be nonnegative at least.
  positivity_indicator <- ifelse(is_positive_component, 0, Inf)
  negativity_indicator <- ifelse(is_negative_component, 0, Inf)

  prob$bx <- rbind(
    blx = c(rep(-Inf, n + 1L), rep(0, 2L * p)),
    bux = c(rep(Inf, n + 1L), negativity_indicator, positivity_indicator)
  )

  # Define the constraint $u = y - M x$
  equality_A <- Matrix::Matrix(cbind(
    rep(0, n),
    diag(n),  # $u$
    M,  # $M x_+$
    -M  # $-M x_-$
  ))
  equality_blc <- y
  equality_buc <- y

  # Define the LASSO constraint $\sum_i |x_i| \le V$
  constraint_indicator <- as.numeric(is_sum_constrained_component)
  lasso_A <- c(rep(0, n + 1L), rep(constraint_indicator, 2L))
  lasso_blc <- -Inf
  lasso_buc <- V

  # Enter the linear constraints
  prob$A <- Matrix::Matrix(rbind(equality_A, lasso_A))
  prob$bc <- rbind(
    blc = c(equality_blc, lasso_blc),
    buc = c(equality_buc, lasso_buc)
  )

  # Solve the problem
  r <- Rmosek::mosek(prob, list(verbose = 0L))

  # Return the solution
  stopifnot(identical(r$response$code, 0))
  xx <- r$sol$itr$xx
  xx[(n + 2L):(n + p + 1L)] - xx[(n + p + 2L):(n + 2L * p + 1L)]
}
