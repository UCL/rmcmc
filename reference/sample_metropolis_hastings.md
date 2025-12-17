# Sample from Metropolis-Hastings kernel.

Sample from Metropolis-Hastings kernel.

## Usage

``` r
sample_metropolis_hastings(
  state,
  target_distribution,
  proposal,
  sample_uniform = stats::runif
)
```

## Arguments

- state:

  Current chain state.

- target_distribution:

  Target stationary distribution for chain. A list with named entries
  `log_density` and `gradient_log_density` corresponding to respectively
  functions for evaluating the logarithm of the (potentially
  unnormalized) density of the target distribution and its gradient. As
  an alternative to `gradient_log_density` an entry
  `value_and_gradient_log_density` may instead be provided which is a
  function returning both the value and gradient of the logarithm of the
  (unnormalized) density of the target distribution as a list under the
  names `value` and `gradient` respectively.

- proposal:

  Proposal distribution object. Must define entries `sample`, a function
  to generate sample from proposal distribution given current chain
  state and `log_density_ratio`, a function to compute log density ratio
  for proposal for a given pair of current and proposed chain states.

- sample_uniform:

  Function which generates a random vector from standard uniform
  distribution given an integer size.

## Value

List with named entries

- `state`: corresponding to new chain state,

- `proposed_state`: corresponding to proposed chain state,

- `statistics`: a list with named entries for statistics of transition,
  here this consisting of a named entry `accept_prob` for the Metropolis
  acceptance probability.
