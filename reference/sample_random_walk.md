# Sample new state from random walk proposal.

Sample new state from random walk proposal.

## Usage

``` r
sample_random_walk(
  state,
  target_distribution,
  scale_and_shape,
  sample_auxiliary = stats::rnorm
)
```

## Arguments

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
