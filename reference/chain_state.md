# Construct a new chain state.

The chain state object provides cached access to target distribution log
density and its gradient.

## Usage

``` r
chain_state(position, momentum = NULL)
```

## Arguments

- position:

  Position component of chain state.

- momentum:

  Momentum component of chain state. Optional.

## Value

New chain state object. A list with entries

- `position`: A zero-argument function to evaluate position vector.

- `momentum`: A zero-argument function to evaluate momentum vector.

- `dimension`: A zero-argument function evaluate dimension of position
  and momentum vectors.

- `update`: A function accepting arguments `position` and `momentum` for
  updating the value of one or both of these state components.

- `copy`: A function for creating a copy of the state object including
  any cached values.

- `log_density`: A function accepting argument `target_distribution` for
  evaluating the log density of the target distribution at the current
  state, with caching of the value to avoid recomputation on subsequent
  calls.

- `gradient_log_density`: A function accepting argument
  `target_distribution` for evaluating the gradient of the log density
  of the target distribution at the current state, with caching of the
  value to avoid recomputation on subsequent calls.

## Examples

``` r
state <- chain_state(c(0.1, -0.5))
target_distribution <- list(
  log_density = function(x) -sum(x^2) / 2,
  gradient_log_density = function(x) -x
)
state$gradient_log_density(target_distribution)
#> [1] -0.1  0.5
```
