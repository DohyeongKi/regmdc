
<!-- README.md is generated from README.Rmd. Please edit that file -->

# regmdc

<!-- badges: start -->
<!-- badges: end -->

regmdc is an R package for fitting nonparametric regression models with
mixed derivative constraints. Available estimation methods are entirely
monotonic regression, Hardy—Krause variation denoising, totally convex
regression, MARS via LASSO, and their generalizations. Entirely
monotonic regression and Hardy—Krause variation denoising are
multivariate generalizations of isotonic regression and total variation
denoising, introduced in Fang et al. (2021). MARS via LASSO is a LASSO
variant of the usual MARS method, introduced in Ki et al. (2021). It can
be thought of not only as a multivariate generalization of locally
adaptive regression splines but also as second-order Hardy—Krause
variation denoising. Totally convex regression is a multivariate
extension of univariate convex regression based on total convexity (see,
e.g., Gal (2008)). A document for totally convex regression is now in
preparation.

regmdc mainly consists of the following two generic functions:

-   `get_lasso_problem_soln()`
-   `compute_fit()`

Given an estimation method, `get_lasso_problem_soln()` finds the
solution to the corresponding constrained LASSO problem. For details on
the LASSO problem of each estimation method, see, for example, Section 3
of Fang et al. (2021) (for entirely monotonic regression and
Hardy—Krause variation denoising) and Section 3 and 6 of Ki et
al. (2021) (for MARS via LASSO). There are some other ongoing research
and working papers, and they would be available in the future.

Given the solution to the LASSO problem obtained from
`get_lasso_problem_soln()`, `compute_fit()` computes the fitted values
of an estimation method.

## Installation

You first need to install
[MOSEK](https://docs.mosek.com/latest/install/installation.html), a
software package for optimization, and then
[Rmosek](https://docs.mosek.com/latest/rmosek/install-interface.html),
the MOSEK interface in R.

Once they are properly installed, you can install our package by running
in R the following commands:

``` r
# install.packages("devtools")
devtools::install_github("DohyeongKi/regmdc")
```

## Examples

Here are some examples showing how regmdc can be used. Please refer to
the document of each function in the package for more details.

``` r
################################# 
# Entirely monotonic regression #
#################################
library(regmdc)

fstar <- function(x) {x[1] + x[2]}  # the underlying function

X_design <- expand.grid(rep(list(seq(0, 9.0/10, length.out = 10L)), 2L))  # a design matrix
# Compute the values of f* at the design points
theta <- apply(X_design, MARGIN = 1L, FUN = fstar)

sigma <- 1  # standard Gaussian noise
y <- theta + sigma * rnorm(nrow(X_design))  # an observation vector

# Obtain the solution to the corresponding constrained LASSO problem
raw_soln <- get_lasso_problem_soln(X_design, y, s = 1L, method = "em")$solution

# Remove the zero components from the solution
threshold <- 1e-4  # a threshold to determine whether each component of the solution is zero
is_nonzero <- (abs(raw_soln) > threshold)
processed_soln <- as.numeric(is_nonzero) * raw_soln
compressed_soln <- processed_soln[processed_soln != 0]

X_eval <- expand.grid(rep(list(seq(0, 1, length.out = 51L)), 2L))  # an evaluation matrix
# Compute the fitted values at the evaluation points
compute_fit(X_eval, X_design, s = 1L, method = "em", compressed_soln, is_nonzero)
```

``` r
#################################### 
# Hardy—Krause variation denoising #
####################################
library(regmdc)

fstar <- function(x) {x[1] - x[2]}  # the underlying function

X_design <- expand.grid(rep(list(seq(0, 9.0/10, length.out = 10L)), 2L))  # a design matrix
# Compute the values of f* at the design points
theta <- apply(X_design, MARGIN = 1L, FUN = fstar)

sigma <- 1  # standard Gaussian noise
y <- theta + sigma * rnorm(nrow(X_design))  # an observation vector

# Obtain the solution to the corresponding constrained LASSO problem
raw_soln <- get_lasso_problem_soln(X_design, y, s = 2L, method = "hk", V = 2)$solution

# Remove the zero components from the solution
threshold <- 1e-4  # a threshold to determine whether each component of the solution is zero
is_nonzero <- (abs(raw_soln) > threshold)
processed_soln <- as.numeric(is_nonzero) * raw_soln
compressed_soln <- processed_soln[processed_soln != 0]

X_eval <- expand.grid(rep(list(seq(0, 1, length.out = 51L)), 2L))  # an evaluation matrix
# Compute the fitted values at the evaluation points
compute_fit(X_eval, X_design, s = 2L, method = "hk", compressed_soln, is_nonzero)
```

``` r
######################################################################################## 
# Generalization of entirely monotonic regression and Hardy—Krause variation denoising #
########################################################################################
library(regmdc)

fstar <- function(x) {x[1] - x[2] + x[1] * x[2]}  # the underlying function

X_design <- expand.grid(rep(list(seq(0, 9.0/10, length.out = 10L)), 2L))  # a design matrix
# Compute the values of f* at the design points
theta <- apply(X_design, MARGIN = 1L, FUN = fstar)

sigma <- 1  # standard Gaussian noise
y <- theta + sigma * rnorm(nrow(X_design))  # an observation vector

# Obtain the solution to the corresponding constrained LASSO problem
raw_soln <- get_lasso_problem_soln(X_design, y, s = 2L, method = "emhk", V = 1,
                                   constrained_interactions = c('1-2'),
                                   positive_interactions = c('1'),
                                   negative_interactions = c('2'))$solution

# Remove the zero components from the solution
threshold <- 1e-4  # a threshold to determine whether each component of the solution is zero
is_nonzero <- (abs(raw_soln) > threshold)
processed_soln <- as.numeric(is_nonzero) * raw_soln
compressed_soln <- processed_soln[processed_soln != 0]

X_eval <- expand.grid(rep(list(seq(0, 1, length.out = 51L)), 2L))  # an evaluation matrix
# Compute the fitted values at the evaluation points
compute_fit(X_eval, X_design, s = 2L, method = "emhk", compressed_soln, is_nonzero)
```

## References

\[1\] Ki, D., Fang, B., and Guntuboyina, A. (2021). MARS via LASSO.
<https://arxiv.org/abs/2111.11694>.

\[2\] Fang, B., Guntuboyina, A., and Sen, B. (2021). Multivariate
extensions of isotonic regression and total variation denoising via
entire monotonicity and Hardy—Krause variation. *The Annals of
Statistics*, **49**(2), 769-792.

\[3\] Gal, S. G. (2008). *Shape-Preserving Approximation by Real and
Complex Polynomials*. Birkhäuser, Boston.
