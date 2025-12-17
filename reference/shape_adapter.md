# Create object to adapt proposal shape.

Create object to adapt proposal shape.

## Usage

``` r
shape_adapter(type = "covariance", kappa = 1)
```

## Arguments

- type:

  Type of shape adapter to use. One of:

  - "variance": Diagonal shape matrix adaptation based on estimates of
    target distribution variances (see
    [`variance_shape_adapter()`](http://github-pages.ucl.ac.uk/rmcmc/reference/variance_shape_adapter.md)),

  - "covariance": Dense shape matrix adaptation based on estimates of
    target distribution covariance matrix (see
    [`covariance_shape_adapter()`](http://github-pages.ucl.ac.uk/rmcmc/reference/covariance_shape_adapter.md)).

- kappa:

  Decay rate exponent in `[0.5, 1]` for adaptation learning rate. Value
  of 1 (default) corresponds to computing empirical (co)variances.

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

## See also

[`variance_shape_adapter()`](http://github-pages.ucl.ac.uk/rmcmc/reference/variance_shape_adapter.md),
[`covariance_shape_adapter()`](http://github-pages.ucl.ac.uk/rmcmc/reference/covariance_shape_adapter.md)

## Examples

``` r
proposal <- barker_proposal()
adapter <- shape_adapter()
adapter$initialize(proposal, chain_state(c(0, 0)))
```
