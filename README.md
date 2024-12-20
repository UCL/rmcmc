
<!-- README.md is generated from README.Rmd. Please edit that file -->

# rmcmc: Robust Markov chain Monte Carlo methods

<!-- badges: start -->

[![R-CMD-check](https://github.com/UCL/rmcmc/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/UCL/rmcmc/actions/workflows/R-CMD-check.yaml)
[![pkgdown](https://github.com/UCL/rmcmc/actions/workflows/pkgdown.yaml/badge.svg)](https://github.com/UCL/rmcmc/actions/workflows/pkgdown.yaml)
[![lint](https://github.com/UCL/rmcmc/actions/workflows/lint.yaml/badge.svg)](https://github.com/UCL/rmcmc/actions/workflows/lint.yaml)
[![pre-commit](https://github.com/UCL/rmcmc/actions/workflows/pre-commit.yaml/badge.svg)](https://github.com/UCL/rmcmc/actions/workflows/pre-commit.yaml)
[![codecov](https://codecov.io/github/UCL/rmcmc/graph/badge.svg?token=PL8557fpgT)](https://codecov.io/github/UCL/rmcmc)
<!-- badges: end -->

`rmcmc` is an R package for simulating Markov chains using the Barker
proposal to compute *Markov chain Monte Carlo* (MCMC) estimates of
expectations with respect to a target distribution on a real-valued
vector space. The Barker proposal, described in Livingstone and Zanella
(2022) <https://doi.org/10.1111/rssb.12482>, is a gradient-based MCMC
algorithm inspired by the Barker accept-reject rule. It combines the
robustness of simpler MCMC schemes such as random-walk Metropolis with
the efficiency of gradient-based algorithms such as Metropolis adjusted
Langevin algorithm.

## Installation

You can install the development version of `rmcmc` like so:

``` r
# install.packages("devtools")
devtools::install_github("UCL/rmcmc")
```

## Examples

The snippet belows shows a basic example of using the package to
generate samples from a normal target distribution with random scales.
Adapters are used to tune the proposal scale to achieve a target average
acceptance probability, and to tune the proposal shape with
per-dimension scale factors based on online estimates of the target
distribution variances.

``` r
library(rmcmc)

set.seed(876287L)
dimension <- 3
scales <- exp(rnorm(dimension))
target_distribution <- list(
  log_density = function(x) -sum((x / scales)^2) / 2,
  gradient_log_density = function(x) -x / scales^2
)
proposal <- barker_proposal()
results <- sample_chain(
  target_distribution = target_distribution,
  initial_state = rnorm(dimension),
  n_warm_up_iteration = 10000,
  n_main_iteration = 10000,
  proposal = proposal,
  adapters = list(scale_adapter(), shape_adapter("variance"))
)
mean_accept_prob <- mean(results$statistics[, "accept_prob"])
adapted_shape <- proposal$parameters()$shape
cat(
  sprintf("Average acceptance probability is %.2f", mean_accept_prob),
  sprintf("True target scales: %s", toString(scales)),
  sprintf("Adapter scale est.: %s", toString(adapted_shape)),
  sep = "\n"
)
#> Average acceptance probability is 0.58
#> True target scales: 1.50538046096953, 1.37774732725824, 0.277038897322645
#> Adapter scale est.: 1.5328097767097, 1.42342707172926, 0.280359693392091
```

As a second example, the snippet below demonstrates sampling from a
two-dimensional banana shaped distribution based on the [Rosenbrock
function](https://en.wikipedia.org/wiki/Rosenbrock_function) and
plotting the generated chain samples. Here we use the default values of
the `proposal` and `adapters` arguments to `sample_chain`, corresponding
respectively to the Barker proposal, and adapters for tuning the
proposal scale to coerce the average acceptance rate using a
dual-averaging algorithm, and for tuning the proposal shape based on an
estimate of the target distribution covariance matrix.

``` r
library(rmcmc)

set.seed(651239L)
target_distribution <- list(
  log_density = function(x) -sum(x^2) / 8 - (x[1]^2 - x[2])^2 - (x[1] - 1)^2 / 10,
  gradient_log_density = function(x) {
    c(
      -x[1] / 4 + 4 * x[1] * (x[2] - x[1]^2) - 0.2 * x[1] + 0.2,
      -x[2] / 4 + 2 * x[1]^2 - 2 * x[2]
    )
  }
)
results <- sample_chain(
  target_distribution = target_distribution,
  initial_state = rnorm(2),
  n_warm_up_iteration = 10000,
  n_main_iteration = 10000,
)
plot(
  results$traces[, "position1"],
  results$traces[, "position2"],
  xlab = expression(x[1]),
  ylab = expression(x[2]),
  col = "#1f77b4",
  pch = 20
)
```

<img src="man/figures/README-banana-samples-1.png" width="100%" />
