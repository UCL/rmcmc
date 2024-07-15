
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

This is a basic example which shows you how to sample from a standard
normal distribution, and compute average acceptance probability over
chain.

``` r
library(rmcmc)

target_distribution <- list(
  log_density = function(x) -sum(x^2) / 2,
  grad_log_density = function(x) -x
)
proposal <- barker_proposal(target_distribution, step_size = 2.5)
n_sample <- 1000
dimension <- 2
set.seed(876287L)
state <- chain_state(rnorm(dimension))
sum_accept_prob <- 0
for (s in 2:n_sample) {
  state_and_statistics <- sample_metropolis_hastings(
    state, target_distribution, proposal
  )
  state <- state_and_statistics$state
  accept_prob <- state_and_statistics$statistics$accept_prob
  sum_accept_prob <- sum_accept_prob + accept_prob
}
mean_accept_prob <- sum_accept_prob / n_sample
message(sprintf("Average acceptance probability is %.2f", mean_accept_prob))
#> Average acceptance probability is 0.40
```
