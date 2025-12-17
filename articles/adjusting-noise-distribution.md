# Adjusting the noise distribution in the Barker proposal

The Barker proposal implementation in the `rmcmc` package in
[`barker_proposal()`](http://github-pages.ucl.ac.uk/rmcmc/reference/barker_proposal.md)
by default uses a standard normal distribution for the auxiliary noise
variables generated in the proposal, as suggested by Livingstone and
Zanella ([2022](#ref-livingstone2022barker)). The analysis in Vogrinc,
Livingstone, and Zanella ([2023](#ref-vogrinc2023optimal)) however
implies that alternative choices of noise distribution can give improved
performance in some situations. This vignette demonstrates how to adjust
the Barker proposal noise distribution.

``` r
library(rmcmc)
```

## Example target distribution

As a simple example of a target distribution, we consider a
$D$-dimensional distribution with product form (independent dimensions)
and ‘hyperbolic’ marginals

$$\pi(x) \propto \exp\left( - \sum\limits_{d = 1}^{D}\left( \delta^{2} + x_{d}^{2} \right)^{\frac{1}{2}} \right).$$
The gradient of the corresponding log density function can be derived as
$$\nabla\left( \log\pi \right)(x)_{d} = - x_{d}/\left( \delta^{2} + x_{d}^{2} \right)^{\frac{1}{2}}.$$
This can be implemented in R as

``` r
delta_sq <- 0.1
target_distribution <- list(
  log_density = function(x) -sum(sqrt(delta_sq + x^2)),
  gradient_log_density = function(x) -x / sqrt(delta_sq + x^2)
)
```

Here we use $\delta^{2} = 0.1$, corresponding to Example 4 in Vogrinc,
Livingstone, and Zanella ([2023](#ref-vogrinc2023optimal)).

## Creating Barker proposal with a custom noise distribution

The
[`barker_proposal()`](http://github-pages.ucl.ac.uk/rmcmc/reference/barker_proposal.md)
function accepts an optional argument `sample_auxiliary`, which can be
used to specify the distribution the auxiliary noise variables are drawn
from. The argument should be passed a function accepting a single
integer argument, corresponding to the target distribution dimension,
and returning a vector of auxiliary noise variable samples, one per
target distribution dimension. By default this function is set to
[`stats::rnorm()`](https://rdrr.io/r/stats/Normal.html), thus using
independent standard normal variates for the auxiliary noise variables.

Below we create a function `sample_bimodal` which instead samples from a
bimodal normal mixture, with two equally-weighted normal components with
means $\pm \left( 1 - \sigma^{2} \right)$ and common variance
$\sigma^{2} = 0.01$.

``` r
sample_bimodal <- function(dimension) {
  sigma <- 0.1
  (
    sample(c(-1, 1), dimension, TRUE) * sqrt(1 - sigma^2) +
      stats::rnorm(dimension) * sigma
  )
}
```

This choice of bimodal distribution for the auxiliary noise variables is
motivated by the analysis in Vogrinc, Livingstone, and Zanella
([2023](#ref-vogrinc2023optimal)), which shows that for target
distributions of product form and using a Barker proposal, the
asymptotic *expect squared jump distance* (ESJD) a measure of chain
sampling performance, is maximised when the auxiliary noise distribution
is chosen to be the Rademacher distribution (the distribution with half
of its mass in atoms at -1 and +1). The Rademacher distribution itself
does not give a practical sampling algorithm, as the resulting Markov
kernel will not be irreducible. As a more practical alternative Vogrinc,
Livingstone, and Zanella ([2023](#ref-vogrinc2023optimal)) therefore
suggests using the above normal mixture bimodal distribution, which can
be considered a relaxation of the Rademacher distribution.

We now create two instances of the the Barker proposal, the first using
the default of a standard normal distribution for the auxiliary noise
variables, and the second using the bimodal distribution.

``` r
normal_noise_proposal <- barker_proposal()
bimodal_noise_proposal <- barker_proposal(sample_auxiliary = sample_bimodal)
```

A convenience function
[`bimodal_barker_proposal()`](http://github-pages.ucl.ac.uk/rmcmc/reference/bimodal_barker_proposal.md)
is also provided in the package to simplify creating a proposal with
bimodal noise distribution as above, with optional `sigma` argument
specifying the standard deviation of the normal components. That is the
second line above can equivalently be written
`bimodal_noise_proposal <- bimodal_barker_proposal()` without the need
to define the `sample_bimodal` function.

## Sampling chains with proposals

Now that we have the target distribution and two proposals specified, we
can sample chains to compare their performance.

``` r
set.seed(7861932L)
dimension <- 100
initial_state <- rnorm(dimension)
adapters <- list(scale_adapter())
n_warm_up_iteration <- 10000
n_main_iteration <- 10000
```

Above we specify some common settings for the
[`sample_chain()`](http://github-pages.ucl.ac.uk/rmcmc/reference/sample_chain.md)
calls, namely the dimension of the target distribution, the initial
chain state (shared for both chains), the adapters to use during the
warm-up iterations (here we adapt only the proposal scale to achieve a
target acceptance rate of 0.57), and the number of warm-up and main
chain iterations.

We can now sample a chain using the bimodal noise proposal variant

``` r
bimodal_noise_results <- sample_chain(
  target_distribution = target_distribution,
  proposal = bimodal_noise_proposal,
  initial_state = initial_state,
  n_warm_up_iteration = n_warm_up_iteration,
  n_main_iteration = n_main_iteration,
  adapters = adapters
)
```

and similarly using the default normal noise proposal

``` r
normal_noise_results <- sample_chain(
  target_distribution = target_distribution,
  proposal = normal_noise_proposal,
  initial_state = initial_state,
  n_warm_up_iteration = n_warm_up_iteration,
  n_main_iteration = n_main_iteration,
  adapters = adapters
)
```

## Comparing performance of proposals

For a realisation of a Markov chain $X_{1:N}$, the ESJD metric used to
define the notion of sampling efficiency optimized over in Vogrinc,
Livingstone, and Zanella ([2023](#ref-vogrinc2023optimal)), can be
estimated as

$$\widehat{\textsf{𝖤𝖲𝖩𝖣}}\left( X_{1:N} \right) = \frac{1}{N - 1}\sum\limits_{n = 1}^{N - 1}\| X_{n + 1} - X_{n}\|^{2}$$
with corresponding R implementation as below

``` r
expected_square_jumping_distance <- function(traces) {
  n_iteration <- nrow(traces)
  mean(rowSums(traces[2:n_iteration, ] - traces[1:n_iteration - 1, ])^2)
}
```

Applying this function to the chain traces generated using the Barker
proposal with bimodal noise,

``` r
cat(
  sprintf(
    "Expected square jumping distance using normal noise is %.2f",
    expected_square_jumping_distance(bimodal_noise_results$traces[, 1:dimension])
  )
)
#> Expected square jumping distance using normal noise is 33.75
```

and Barker proposal with normal noise,

``` r
cat(
  sprintf(
    "Expected square jumping distance using bimodal noise is %.2f",
    expected_square_jumping_distance(normal_noise_results$traces[, 1:dimension])
  )
)
#> Expected square jumping distance using bimodal noise is 17.42
```

in both cases selecting columns `1:dimension` to exclude the traced log
density values, we find that the proposal with bimodal noise gives
roughly a factor two improvement in ESJD compared to using normal noise,
correlating with the results in Vogrinc, Livingstone, and Zanella
([2023](#ref-vogrinc2023optimal)).

``` r
library(posterior)
#> This is posterior version 1.6.1
#> 
#> Attaching package: 'posterior'
#> The following objects are masked from 'package:stats':
#> 
#>     mad, sd, var
#> The following objects are masked from 'package:base':
#> 
#>     %in%, match
```

Using the `posterior` R package, we can all compare how the chains
perform in terms of the estimated *effective sample size* (ESS) of the
estimate of the mean of the target log density.

First considering the chain using the Barker proposal with bimodal noise

``` r
cat(
  sprintf(
    "Estimated ESS of mean(target_log_density) using bimodal noise is %.0f",
    ess_mean(
      extract_variable(
        as_draws_matrix(bimodal_noise_results$traces), "target_log_density"
      )
    )
  )
)
#> Estimated ESS of mean(target_log_density) using bimodal noise is 345
```

and similarly for the Barker proposal with normal noise

``` r
cat(
  sprintf(
    "Estimated ESS of mean(target_log_density) using normal noise is %.0f",
    ess_mean(
      extract_variable(
        as_draws_matrix(normal_noise_results$traces), "target_log_density"
      )
    )
  )
)
#> Estimated ESS of mean(target_log_density) using normal noise is 274
```

we again find that the bimodal noise variant appears to show
significantly improved sampling efficiency.

## References

Livingstone, Samuel, and Giacomo Zanella. 2022. “The Barker proposal:
combining robustness and efficiency in gradient-based MCMC.” *Journal of
the Royal Statistical Society Series B: Statistical Methodology* 84 (2):
496–523. <https://doi.org/10.1111/rssb.12482>.

Vogrinc, Jure, Samuel Livingstone, and Giacomo Zanella. 2023. “Optimal
Design of the Barker Proposal and Other Locally Balanced
Metropolis–Hastings Algorithms.” *Biometrika* 110 (3): 579–95.
