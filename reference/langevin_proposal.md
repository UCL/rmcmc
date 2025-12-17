# Create a new Langevin proposal object.

The Langevin proposal is a gradient-based proposal corresponding to a
Euler-Maruyama time discretisation of a Langevin diffusion.

## Usage

``` r
langevin_proposal(scale = NULL, shape = NULL, sample_auxiliary = stats::rnorm)
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

## References

Besag, J. (1994). "Comments on "Representations of knowledge in complex
systems" by U. Grenander and MI Miller". *Journal of the Royal
Statistical Society, Series B*. 56: 591–592.

Roberts, G. O., & Tweedie, R. L. (1996). Exponential convergence of
Langevin distributions and their discrete approximations. *Bernoulli* 2
(4), 341 - 363.

## Examples

``` r
target_distribution <- list(
  log_density = function(x) -sum(x^2) / 2,
  gradient_log_density = function(x) -x
)
proposal <- langevin_proposal(scale = 1.)
state <- chain_state(c(0., 0.))
withr::with_seed(
  876287L, proposed_state <- proposal$sample(state, target_distribution)
)
log_density_ratio <- proposal$log_density_ratio(
  state, proposed_state, target_distribution
)
proposal$update(scale = 0.5)
```
