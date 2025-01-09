
<!-- README.md is generated from README.Rmd. Please edit that file -->

# regmdc

<!-- badges: start -->
<!-- badges: end -->

regmdc is an R package for fitting nonparametric regression models with
mixed derivative constraints. Available estimation methods are entirely
monotonic regression, Hardy—Krause variation denoising, totally concave
regression, MARS via LASSO, and their generalizations. Entirely
monotonic regression and Hardy—Krause variation denoising are
multivariate generalizations of isotonic regression and total variation
denoising, introduced in Fang et al. (2021). MARS via LASSO is a LASSO
variant of the usual MARS method, introduced in Ki et al. (2024). It can
be thought of not only as a multivariate generalization of locally
adaptive regression splines but also as second-order Hardy—Krause
variation denoising. Totally concave regression is a multivariate
extension of univariate concave regression based on total concavity
(see, e.g., Gal (2008)), proposed in Ki and Guntuboyina (2025+).

regmdc mainly consists of the following two generic functions:

- `regmdc()`
- `predict_regmdc()`

Given an estimation method, `regmdc()` builds the model fit to data, by
solving the corresponding constrained LASSO problem. For details on the
LASSO problem of each estimation method, see, for example, Section 3 of
Fang et al. (2021) (for entirely monotonic regression and Hardy—Krause
variation denoising), Section 2 of Ki et al. (2024) (for MARS via
LASSO), and Section 3 of Ki and Guntuboyina (2025+) (for totally concave
regression).

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
library(regmdc)
################################################################################ 
# (1) Entirely monotonic regression    
# (2) Hardy—Krause variation denoising 
# (3) Generalization of entirely monotonic regression and Hardy—Krause variation denoising 
################################################################################
fstar <- function(x) {(
  (x[1] - 0.25 >= 0) + (x[2] - 0.25 >= 0) 
  + (x[1] - 0.25 >= 0) * (x[2] - 0.25 >= 0)
)}  # the true underlying function
X_design <- expand.grid(rep(list(seq(0, 1, length.out = 10L)), 3L))  # a design matrix
colnames(X_design) <- c("VarA", "VarB", "VarC")
theta <- apply(X_design, MARGIN = 1L, FUN = fstar)  # the values of f* at the design points
sigma <- 0.1  # standard Gaussian noise
y <- theta + sigma * rnorm(nrow(X_design))  # an observation vector of a response variable

# Build an entirely monotonic regression model
# em_model <- regmdc(X_design, y, s = 2L, method = "em")
# em_model <- regmdc(X_design, y, s = 2L, method = "em", is_scaled = TRUE)
em_model <- regmdc(X_design, y, s = 2L, method = "em", is_lattice = TRUE)
# em_model <- regmdc(X_design, y, s = 2L, method = "em",
#                    is_monotonically_increasing = FALSE)
# em_model <- regmdc(X_design, y, s = 2L, method = "em",
#                    increasing_covariates = c(1L, 2L))
# em_model <- regmdc(X_design, y, s = 2L, method = "em",
#                    decreasing_covariates = c(3L))
# em_model <- regmdc(X_design, y, s = 2L, method = "em",
#                    increasing_covariates = c("VarA", "VarB"),
#                    decreasing_covariates = c("VarC"))

# Build a Hardy-Krause variation denoising model
# hk_model <- regmdc(X_design, y, s = 2L, method = "hk", V = 3.0)
# hk_model <- regmdc(X_design, y, s = 2L, method = "hk", V = 3.0, is_scaled = TRUE)
hk_model <- regmdc(X_design, y, s = 2L, method = "hk", V = 3.0, is_lattice = TRUE)

# Build a generalized model
# emhk_model <- regmdc(X_design, y, s = 2L, method = "emhk", V = 2.0, is_lattice = TRUE,
#                      variation_constrained_covariates = c(2L))
# emhk_model <- regmdc(X_design, y, s = 2L, method = "emhk", V = 2.0,
#                      is_monotonically_increasing = FALSE,
#                      variation_constrained_covariates = c("VarB"))
emhk_model <- regmdc(X_design, y, s = 2L, method = "emhk", V = 2.0,
                     increasing_covariates = c(1L),
                     variation_constrained_covariates = c(2L))
# emhk_model <- regmdc(X_design, y, s = 2L, method = "emhk", V = 2.0,
#                      increasing_covariates = c("VarA"),
#                      decreasing_covariates = c("VarC"),
#                      variation_constrained_covariates = c("VarB"))

# Generate predictions at new data points
X_pred <- c(1.0/3, 2.0/3, 1.0/3)
predict_regmdc(em_model, X_pred)
predict_regmdc(hk_model, X_pred)
predict_regmdc(emhk_model, X_pred)
X_pred <- matrix(c(1.0/3, 2.0/3, 1.0/3, 
                   2.0/3, 1.0/3, 2.0/3), 
                 ncol = 3L, byrow = TRUE)
predict_regmdc(em_model, X_pred)
predict_regmdc(hk_model, X_pred)
predict_regmdc(emhk_model, X_pred)
```

``` r
library(regmdc)
################################################################################ 
# (4) Totally concave regression  
# (5) MARS via LASSO
# (6) Generalization of totally concave regression and MARS via LASSO 
################################################################################
fstar <- function(x) {(
  - max(x[1] - 0.25, 0) - max(x[2] - 0.25, 0)
  - max(x[1] - 0.25, 0) * max(x[2] - 0.25, 0)
)}  # the true underlying function

X_design <- cbind(runif(100), runif(100), runif(100))
colnames(X_design) <- c("VarA", "VarB", "VarC")
theta <- apply(X_design, MARGIN = 1L, FUN = fstar)  # the values of f* at the design points
sigma <- 0.1  # standard Gaussian noise
y <- theta + sigma * rnorm(nrow(X_design))  # an observation vector

# Build a totally convex regression model
tc_model <- regmdc(X_design, y, s = 2L, method = "tc")
# tc_model <- regmdc(X_design, y, s = 2L, method = "tc",
#                    is_totally_concave = FALSE)
# tc_model <- regmdc(X_design, y, s = 2L, method = "tc",
#                    concave_covariates = c(1L, 2L))
# tc_model <- regmdc(X_design, y, s = 2L, method = "tc",
#                    convex_covariates = c(3L))
# tc_model <- regmdc(X_design, y, s = 2L, method = "tc",
#                    concave_covariates = c("VarA", "VarB"),
#                    convex_covariates = c("VarC"))
# tc_model <- regmdc(X_design, y, s = 2L, method = "tc",
#                    extra_linear_covariates = c(3L))
# tc_model <- regmdc(X_design, y, s = 2L, method = "tc",
#                    is_totally_concave = FALSE,
#                    extra_linear_covariates = c("VarC"))
# tc_model <- regmdc(X_design, y, s = 2L, method = "tc",
#                    extra_linear_covariates = c(2L, 3L))
# tc_model <- regmdc(X_design, y, s = 2L, method = "tc",
#                    concave_covariates = c("VarA"),
#                    extra_linear_covariates = c("VarC"))
# tc_model <- regmdc(X_design, y, s = 2L, method = "tc",
#                    number_of_bins = 20L,
#                    extra_linear_covariates = c(3L))

# Build a MARS via LASSO model
mars_model <- regmdc(X_design, y, s = 2L, method = "mars", V = 3.0)
# mars_model <- regmdc(X_design, y, s = 2L, method = "mars", V = 3.0,
#                      number_of_bins = 20L)
# mars_model <- regmdc(X_design, y, s = 2L, method = "mars", V = 3.0,
#                      number_of_bins = c(10L, 20L, 20L))
# mars_model <- regmdc(X_design, y, s = 2L, method = "mars", V = 3.0,
#                      number_of_bins = c(10L, 20L, NA))
# mars_model <- regmdc(X_design, y, s = 2L, method = "mars", V = 3.0,
#                      number_of_bins = c(10L, 20L, NA),
#                      extra_linear_covariates = c("VarC"))

# Build a generalized model
# tcmars_model <- regmdc(X_design, y, s = 2L, method = "tcmars", V = 2.0,
#                        variation_constrained_covariates = c(2L))
# tcmars_model <- regmdc(X_design, y, s = 2L, method = "tcmars", V = 2.0,
#                        is_totally_concave = FALSE,
#                        variation_constrained_covariates = c(2L))
tcmars_model <- regmdc(X_design, y, s = 2L, method = "tcmars", V = 2.0,
                       concave_covariates = c(1L),
                       variation_constrained_covariates = c(2L))
# tcmars_model <- regmdc(X_design, y, s = 2L, method = "tcmars", V = 2.0,
#                        concave_covariates = c("VarA"),
#                        convex_covariates = c("VarC"),
#                        variation_constrained_covariates = c("VarB"))
# tcmars_model <- regmdc(X_design, y, s = 2L, method = "tcmars", V = 2.0,
#                        concave_covariates = c(1L),
#                        variation_constrained_covariates = c(2L),
#                        extra_linear_covariates = c(3L))
# tcmars_model <- regmdc(X_design, y, s = 2L, method = "tcmars", V = 2.0, 
#                        number_of_bins = 20L,
#                        concave_covariates = c("VarA"),
#                        variation_constrained_covariates = c("VarB"),
#                        extra_linear_covariates = c("VarC"))

# Generate predictions at new data points
X_pred <- c(1.0/3, 2.0/3, 1.0/3)
predict_regmdc(tc_model, X_pred)
predict_regmdc(mars_model, X_pred)
predict_regmdc(tcmars_model, X_pred)
X_pred <- matrix(c(1.0/3, 2.0/3, 1.0/3, 
                   2.0/3, 1.0/3, 2.0/3), 
                 ncol = 3L, byrow = TRUE)
predict_regmdc(tc_model, X_pred)
predict_regmdc(mars_model, X_pred)
predict_regmdc(tcmars_model, X_pred)
```

## References

\[1\] Ki, D. and Guntuboyina, A. (2025+). Totally Concave Regression.
Available at <https://arxiv.org/abs/2501.04360>.

\[2\] Ki, D., Fang, B., and Guntuboyina, A. (2024). MARS via LASSO.
*Annals of Statistics*, **52**(3), 1102-1126.

\[3\] Fang, B., Guntuboyina, A., and Sen, B. (2021). Multivariate
extensions of isotonic regression and total variation denoising via
entire monotonicity and Hardy—Krause variation. *Annals of Statistics*,
**49**(2), 769-792.

\[4\] Gal, S. G. (2008). *Shape-Preserving Approximation by Real and
Complex Polynomials*. Birkhäuser, Boston.
