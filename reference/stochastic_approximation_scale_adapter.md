# Create object to adapt proposal scale to coerce average acceptance rate using a Robbins and Monro (1951) scheme.

When combined with
[`covariance_shape_adapter()`](http://github-pages.ucl.ac.uk/rmcmc/reference/covariance_shape_adapter.md)
corresponds to Algorithm 4 in Andrieu and Thoms (2009).

## Usage

``` r
stochastic_approximation_scale_adapter(
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

Andrieu, C., & Thoms, J. (2008). A tutorial on adaptive MCMC.
*Statistics and Computing*, 18, 343-373.

Robbins, H., & Monro, S. (1951). A stochastic approximation method. *The
Annals of Mathematical Statistics*, 400-407.

## Examples

``` r
proposal <- barker_proposal()
adapter <- stochastic_approximation_scale_adapter(
  initial_scale = 1., target_accept_prob = 0.4
)
adapter$initialize(proposal, chain_state(c(0, 0)))
```
