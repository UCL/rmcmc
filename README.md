
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

## Example

This is a basic example which shows you how to generate samples from a
normal target distribution with random scales. Adapters are used to tune
the proposal scale to achieve a target average acceptance probability;
to tune the proposal shape with per-dimension scale factors based on
online estimates of the target distribution variances.

``` r
library(rmcmc)

set.seed(876287L)
dimension <- 3
scales <- exp(rnorm(dimension))
target_distribution <- list(
  log_density = function(x) -sum((x / scales)^2) / 2,
  gradient_log_density = function(x) -x / scales^2
)
proposal <- barker_proposal(target_distribution)
results <- sample_chain(
  target_distribution = target_distribution,
  proposal = proposal,
  initial_state = rnorm(dimension),
  n_warm_up_iteration = 1000,
  n_main_iteration = 1000,
  adapters = list(simple_scale_adapter(), variance_shape_adapter())
)
mean_accept_prob <- mean(results$statistics[, "accept_prob"])
adapted_shape <- proposal$parameters()$shape
cat(
  sprintf("Average acceptance probability is %.2f", mean_accept_prob),
  sprintf("True target scales: %s", toString(scales)),
  sprintf("Adapter scale est.: %s", toString(adapted_shape)),
  sep = "\n"
)
#> Average acceptance probability is 0.41
#> True target scales: 1.50538046096953, 1.37774732725824, 0.277038897322645
#> Adapter scale est.: 1.2489768457131, 1.23111560302158, 0.215024121396933
```
