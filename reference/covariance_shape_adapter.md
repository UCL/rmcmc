# Create object to adapt proposal with shape based on estimate of target distribution covariance matrix.

Corresponds to Algorithm 2 in Andrieu and Thoms (2009), which is itself
a restatement of method proposed in Haario et al. (2001).

## Usage

``` r
covariance_shape_adapter(kappa = 1)
```

## Arguments

- kappa:

  Decay rate exponent in `[0.5, 1]` for adaptation learning rate. Value
  of 1 (default) corresponds to computing empirical covariance matrix.

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

## Details

Requires `ramcmc` package to be installed for access to efficient rank-1
Cholesky update function
[`ramcmc::chol_update`](https://rdrr.io/pkg/ramcmc/man/chol_update.html).

## References

Andrieu, C., & Thoms, J. (2008). A tutorial on adaptive MCMC.
*Statistics and Computing*, 18, 343-373.

Haario, H., Saksman, E., & Tamminen, J. (2001). An adaptive Metropolis
algorithm. *Bernoulli*, 7(2): 223-242.

## Examples

``` r
proposal <- barker_proposal()
adapter <- covariance_shape_adapter()
adapter$initialize(proposal, chain_state(c(0, 0)))
```
