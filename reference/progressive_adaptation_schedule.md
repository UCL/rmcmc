# Create a progressive adaptation schedule for use with [`sample_chain()`](http://github-pages.ucl.ac.uk/rmcmc/reference/sample_chain.md).

Returns a function (a schedule constructor) that, given the total number
of warm-up iterations, produces a three-stage adaptation schedule:

## Usage

``` r
progressive_adaptation_schedule(
  n_fixed_shape_iteration = 50L,
  n_diagonal_shape_iteration = 50L,
  scale_adapter_ = scale_adapter(),
  diagonal_shape_adapter_ = shape_adapter("variance"),
  dense_shape_adapter_ = shape_adapter("covariance")
)
```

## Arguments

- n_fixed_shape_iteration:

  Number of iterations in Stage 1 (scale only). Default 50.

- n_diagonal_shape_iteration:

  Number of iterations in Stage 2 (scale + diagonal shape). Default 50.

- scale_adapter\_:

  Adapter object for the scale. Defaults to
  [`scale_adapter()`](http://github-pages.ucl.ac.uk/rmcmc/reference/scale_adapter.md).

- diagonal_shape_adapter\_:

  Adapter object for the diagonal shape stage. Defaults to
  `shape_adapter("variance")`.

- dense_shape_adapter\_:

  Adapter object for the dense shape stage. Defaults to
  `shape_adapter("covariance")`.

## Value

A function that accepts `n_warm_up_iteration` and returns a staged
adapter list suitable for the `adapters` argument of
[`sample_chain()`](http://github-pages.ucl.ac.uk/rmcmc/reference/sample_chain.md).

## Details

- **Stage 1** (`n_fixed_shape_iteration` iterations): adapt scale only,
  keeping the proposal shape fixed at the identity. This lets the step
  size stabilise before any covariance estimation begins.

- **Stage 2** (`n_diagonal_shape_iteration` iterations): adapt scale and
  learn a diagonal proposal shape (per-dimension variances).

- **Stage 3** (remaining iterations): adapt scale and learn a dense
  proposal shape (full covariance matrix).

If the sum of `n_fixed_shape_iteration` and `n_diagonal_shape_iteration`
exceeds the total number of warm-up iterations, the two counts are
reduced proportionally so that each gets half of the available
iterations and Stage 3 is skipped.

## Examples

``` r
target_distribution <- list(
  log_density = function(x) -sum(x^2) / 2,
  gradient_log_density = function(x) -x
)
withr::with_seed(876287L, {
  results <- sample_chain(
    target_distribution,
    initial_state = stats::rnorm(2),
    n_warm_up_iteration = 1000,
    n_main_iteration = 1000,
    adapters = progressive_adaptation_schedule()
  )
})
```
