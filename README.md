
<!-- README.md is generated from README.Rmd. Please edit that file -->

# rmcmc: Robust Markov chain Monte Carlo methods

<!-- badges: start -->

[![R-CMD-check](https://github.com/UCL/rmcmc/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/UCL/rmcmc/actions/workflows/R-CMD-check.yaml)
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
dimension <- 2
scales <- exp(2 * rnorm(dimension))
target_distribution <- list(
  log_density = function(x) -sum((x / scales)^2) / 2,
  grad_log_density = function(x) -x / scales^2
)
proposal <- barker_proposal(target_distribution)
adapters <- list(
  scale_adapter(proposal, initial_scale = 1., target_accept_prob = 0.4),
  variance_adapter(proposal)
)
n_sample <- 1000
state <- chain_state(rnorm(dimension))
for (adapter in adapters) {
  adapter$initialize(state)
}
sum_accept_prob <- 0
for (s in 2:n_sample) {
  state_and_statistics <- sample_metropolis_hastings(
    state, target_distribution, proposal
  )
  for (adapter in adapters) {
    adapter$update(s, state_and_statistics)
  }
  state <- state_and_statistics$state
  accept_prob <- state_and_statistics$statistics$accept_prob
  sum_accept_prob <- sum_accept_prob + accept_prob
}
mean_accept_prob <- sum_accept_prob / n_sample
adapted_shape <- proposal$parameters()$shape
cat(
  sprintf("Average acceptance probability is %.2f", mean_accept_prob),
  sprintf("True target scales: %s", toString(scales)),
  sprintf("Adapter scale est.: %s", toString(adapted_shape)),
  sep = "\n"
)
#> Average acceptance probability is 0.39
#> True target scales: 2.26617033226883, 1.89818769776724
#> Adapter scale est.: 2.21142958110019, 1.68681689502187
```
