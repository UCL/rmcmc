# Construct target distribution from a BridgeStan `StanModel` object.

Construct target distribution from a BridgeStan `StanModel` object.

## Usage

``` r
target_distribution_from_stan_model(
  model,
  include_log_density = TRUE,
  include_generated_quantities = FALSE,
  include_transformed_parameters = FALSE,
  seed = 1234L
)
```

## Arguments

- model:

  Stan model object to use for target (posterior) distribution.

- include_log_density:

  Whether to include an entry `log_density` corresponding to current log
  density for target distribution in values returned by trace function.

- include_generated_quantities:

  Whether to included generated quantities in Stan model definition in
  values returned by trace function.

- include_transformed_parameters:

  Whether to include transformed parameters in Stan model definition in
  values returned by trace function.

- seed:

  Seed to use for random number generator used to for any random numbers
  generated as part of generated quantities block.

## Value

A list with entries

- `log_density`: A function to evaluate log density function for target
  distribution given current position vector.

- `value_and_gradient_log_density`: A function to evaluate value and
  gradient of log density function for target distribution given current
  position vector, returning as a list with entries `value` and
  `gradient`.

- `trace_function`: A function which given a
  [`chain_state()`](http://github-pages.ucl.ac.uk/rmcmc/reference/chain_state.md)
  object returns a named vector of values to trace during sampling. The
  constrained parameter values of model will always be included.

## Examples

``` r
model <- example_gaussian_stan_model()
target_distribution <- target_distribution_from_stan_model(model)
withr::with_seed(
  876287L, state <- chain_state(stats::rnorm(model$param_unc_num()))
)
state$log_density(target_distribution)
#> [1] -18.93961
target_distribution$trace_function(state)
#>          mu       sigma log_density 
#>  -0.5220125   0.9487576 -18.9396077 
```
