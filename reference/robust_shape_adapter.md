# Create object to adapt proposal shape (and scale) using robust adaptive Metropolis algorithm of Vihola (2012).

Requires `ramcmc` package to be installed.

## Usage

``` r
robust_shape_adapter(
  initial_scale = NULL,
  target_accept_prob = NULL,
  kappa = 0.6
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

Vihola, M. (2012). Robust adaptive Metropolis algorithm with coerced
acceptance rate. *Statistics and Computing*, 22, 997-1008.
[doi:10.1007/s11222-011-9269-5](https://doi.org/10.1007/s11222-011-9269-5)

## Examples

``` r
proposal <- barker_proposal()
adapter <- robust_shape_adapter(initial_scale = 1., target_accept_prob = 0.4)
adapter$initialize(proposal, chain_state(c(0, 0)))
```
