
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

-   `regmdc()`
-   `predict_regmdc()`

Given an estimation method, `regmdc()` builds the model fit to data, by
solving the corresponding constrained LASSO problem. For details on the
LASSO problem of each estimation method, see, for example, Section 3 of
Fang et al. (2021) (for entirely monotonic regression and Hardy—Krause
variation denoising) and Section 3 and 6 of Ki et al. (2021) (for MARS
via LASSO). There are some other ongoing research and working papers,
and they would be available in the future.

Given the model obtained from `regmdc()`, `predict_regmdc()` provides
predictions at new data points.

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

fstar <- function(x) {x[1] + x[2]}  # the true underlying function
X_design <- expand.grid(rep(list(seq(0, 9.0/10, length.out = 10L)), 2L))  # a design matrix
theta <- apply(X_design, MARGIN = 1L, FUN = fstar)  # the values of f* at the design points
sigma <- 1.0  # standard Gaussian noise
y <- theta + sigma * rnorm(nrow(X_design))  # an observation vector

# Build an entirely monotonic regression model
em_model <- regmdc(X_design, y, s = 1L, method = "em", threshold = 1e-04)

# Generate predictions at new data points
X_pred <- matrix(c(1.0/3, 2.0/3, 2.0/3, 1.0/3), nrow = 2L, ncol = 2L)
predict_regmdc(em_model, X_pred)
```

``` r
#################################### 
# Hardy—Krause variation denoising #
####################################
library(regmdc)

fstar <- function(x) {x[1] - x[2]}  # the true underlying function
X_design <- expand.grid(rep(list(seq(0, 9.0/10, length.out = 10L)), 2L))  # a design matrix
theta <- apply(X_design, MARGIN = 1L, FUN = fstar)  # the values of f* at the design points
sigma <- 1.0  # standard Gaussian noise
y <- theta + sigma * rnorm(nrow(X_design))  # an observation vector

# Build a Hardy-Krause variation denoising model
hk_model <- regmdc(X_design, y, s = 2L, method = "hk", V = 2.0, threshold = 1e-04)

# Generate predictions at new data points
X_pred <- matrix(c(1.0/3, 2.0/3, 2.0/3, 1.0/3), nrow = 2L, ncol = 2L)
predict_regmdc(hk_model, X_pred)
```

``` r
######################################################################################## 
# Generalization of entirely monotonic regression and Hardy—Krause variation denoising #
########################################################################################
library(regmdc)

fstar <- function(x) {x[1] - x[2] + x[1] * x[2]}  # the true underlying function
X_design <- expand.grid(rep(list(seq(0, 9.0/10, length.out = 10L)), 2L))  # a design matrix
theta <- apply(X_design, MARGIN = 1L, FUN = fstar)  # the values of f* at the design points
sigma <- 1.0  # standard Gaussian noise
y <- theta + sigma * rnorm(nrow(X_design))  # an observation vector

# Build a generalized model
emhk_model <- regmdc(X_design, y, s = 2L, method = "emhk", V = 1.0, 
                     threshold = 1e-04, 
                     constrained_interactions = c('1-2'),
                     positive_interactions = c('1'),
                     negative_interactions = c('2'))

# Generate predictions at new data points
X_pred <- matrix(c(1.0/3, 2.0/3, 2.0/3, 1.0/3), nrow = 2L, ncol = 2L)
predict_regmdc(emhk_model, X_pred)
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
