# Create object to adapt proposal scale to coerce average acceptance rate.

Create object to adapt proposal scale to coerce average acceptance rate.

## Usage

``` r
scale_adapter(
  algorithm = "dual_averaging",
  initial_scale = NULL,
  target_accept_prob = NULL,
  ...
)
```

## Arguments

- algorithm:

  String specifying algorithm to use. One of:

  - "stochastic_approximation" to use a Robbins-Monro (1951) based
    scheme,

  - "dual_averaging" to use dual-averaging scheme of Nesterov (2009).

- initial_scale:

  Initial value to use for scale parameter. If not set explicitly a
  proposal and dimension dependent default will be used.

- target_accept_prob:

  Target value for average accept probability for chain. If not set a
  proposal dependent default will be used.

- ...:

  Any additional algorithmic parameters to pass through to
  [`dual_averaging_scale_adapter()`](http://github-pages.ucl.ac.uk/rmcmc/reference/dual_averaging_scale_adapter.md)
  or
  [`stochastic_approximation_scale_adapter()`](http://github-pages.ucl.ac.uk/rmcmc/reference/stochastic_approximation_scale_adapter.md).

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

Robbins, H., & Monro, S. (1951). A stochastic approximation method. *The
Annals of Mathematical Statistics*, 400-407.

## See also

[`dual_averaging_scale_adapter()`](http://github-pages.ucl.ac.uk/rmcmc/reference/dual_averaging_scale_adapter.md),
[`stochastic_approximation_scale_adapter()`](http://github-pages.ucl.ac.uk/rmcmc/reference/stochastic_approximation_scale_adapter.md)

## Examples

``` r
proposal <- barker_proposal()
adapter <- scale_adapter(initial_scale = 1., target_accept_prob = 0.4)
adapter$initialize(proposal, chain_state(c(0, 0)))
```
