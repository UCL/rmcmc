
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

This is a basic example which shows you how to solve a common problem:

``` r
library(rmcmc)
## basic example code
```
