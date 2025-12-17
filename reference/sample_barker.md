# Sample new state from Barker proposal.

Sample new state from Barker proposal.

## Usage

``` r
sample_barker(
  state,
  target_distribution,
  scale_and_shape,
  sample_auxiliary = stats::rnorm,
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

- scale_and_shape:

  Scalar, vector or matrix which scales and shapes proposal
  distribution. If a scalar (in which case the value should be
  non-negative) the auxiliary vector will be isotropically scaled by the
  value. If a vector (in which case the value should be equal in length
  to the dimension of the space and all entries non-negative) each
  dimension of the auxiliary vector will be scaled separately. If a
  matrix (in which case the value should be a square matrix with size
  equal to the dimension of the space) then by pre-multiplying the
  auxiliary vector arbitrary linear transformations can be performed.

- sample_auxiliary:

  Function which generates a random vector from auxiliary variable
  distribution.

- sample_uniform:

  Function which generates a random vector from standard uniform
  distribution given an integer size.

## Value

Proposed new chain state.
