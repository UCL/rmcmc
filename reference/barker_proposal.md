# Create a new Barker proposal object.

The Barker proposal is a gradient-based proposal inspired by the Barker
accept-reject rule and proposed in Livingstone and Zanella (2022). It
offers improved robustness compared to alternative gradient-based
proposals such as Langevin proposals.

## Usage

``` r
barker_proposal(
  scale = NULL,
  shape = NULL,
  sample_auxiliary = stats::rnorm,
  sample_uniform = stats::runif
)
```

## Arguments

- scale:

  Scale parameter of proposal distribution. A non-negative scalar value
  determining scale of steps proposed.

- shape:

  Shape parameter of proposal distribution. Either a vector
  corresponding to a diagonal shape matrix with per-dimension scaling
  factors, or a matrix allowing arbitrary linear transformations.

- sample_auxiliary:

  Function which generates a random vector from auxiliary variable
  distribution.

- sample_uniform:

  Function which generates a random vector from standard uniform
  distribution given an integer size.

## Value

Proposal object. A list with entries

- `sample`: a function to generate sample from proposal distribution
  given current chain state,

- `log_density_ratio`: a function to compute log density ratio for
  proposal for a given pair of current and proposed chain states,

- `update`: a function to update parameters of proposal,

- `parameters`: a function to return list of current parameter values.

- `default_target_accept_prob`: a function returning the default target
  acceptance rate to use for any scale adaptation.

- `default_initial_scale`: a function which given a dimension gives a
  default value to use for the initial proposal scale parameter.

## Details

For more details see the vignette:
[`vignette("barker-proposal", package = "rmcmc")`](http://github-pages.ucl.ac.uk/rmcmc/articles/barker-proposal.md)

## References

Livingstone, S., & Zanella, G. (2022). The Barker proposal: combining
robustness and efficiency in gradient-based MCMC. *Journal of the Royal
Statistical Society Series B: Statistical Methodology*, 84(2), 496-523.

## Examples

``` r
target_distribution <- list(
  log_density = function(x) -sum(x^2) / 2,
  gradient_log_density = function(x) -x
)
proposal <- barker_proposal(scale = 1.)
state <- chain_state(c(0., 0.))
withr::with_seed(
  876287L, proposed_state <- proposal$sample(state, target_distribution)
)
log_density_ratio <- proposal$log_density_ratio(
  state, proposed_state, target_distribution
)
proposal$update(scale = 0.5)
```
