# Create object to adapt proposal scale to coerce average acceptance rate using dual averaging scheme of Nesterov (2009) and Hoffman and Gelman (2014).

Create object to adapt proposal scale to coerce average acceptance rate
using dual averaging scheme of Nesterov (2009) and Hoffman and Gelman
(2014).

## Usage

``` r
dual_averaging_scale_adapter(
  initial_scale = NULL,
  target_accept_prob = NULL,
  kappa = 0.75,
  gamma = 0.05,
  iteration_offset = 10,
  mu = NULL
)
```

## Arguments

- initial_scale:

  Initial value to use for scale parameter. If not set explicitly a
  proposal and dimension dependent default will be used.

- target_accept_prob:

  Target value for average accept probability for chain. If not set a
  proposal dependent default will be used.

- kappa:

  Decay rate exponent in `[0.5, 1]` for adaptation learning rate.
  Defaults to value recommended in Hoffman and Gelman (2014).

- gamma:

  Regularization coefficient for (log) scale in dual averaging
  algorithm. Controls amount of regularization of (log) scale towards
  `mu`. Should be set to a non-negative value. Defaults to value
  recommended in Hoffman and Gelman (2014).

- iteration_offset:

  Offset to chain iteration used for the iteration based weighting of
  the adaptation statistic error estimate. Should be set to a
  non-negative value. A value greater than zero has the effect of
  stabilizing early iterations. Defaults to value recommended in Hoffman
  and Gelman (2014).

- mu:

  Value to regularize (log) scale towards. If `NULL` (the default), `mu`
  will be set to `log(10 * initial_scale)`, as recommended in Hoffman
  and Gelman (2014).

## Value

List of functions with entries

- `initialize`, a function for initializing adapter state and proposal
  parameters at beginning of chain,

- `update` a function for updating adapter state and proposal parameters
  on each chain iteration,

- `finalize` a function for performing any final updates to adapter
  state and proposal parameters on completion of chain sampling (may be
  `NULL` if unused).

- `state` a zero-argument function for accessing current values of
  adapter state variables.

## References

Nesterov, Y. (2009). Primal-dual subgradient methods for convex
problems. *Mathematical Programming*, 120(1), 221-259.

Hoffman, M. D., & Gelman, A. (2014). The No-U-Turn sampler: adaptively
setting path lengths in Hamiltonian Monte Carlo. *Journal of Machine
Learning Research*, 15(1), 1593-1623.

## Examples

``` r
proposal <- barker_proposal()
adapter <- dual_averaging_scale_adapter(
  initial_scale = 1., target_accept_prob = 0.4
)
adapter$initialize(proposal, chain_state(c(0, 0)))
```
