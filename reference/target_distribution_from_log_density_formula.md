# Construct target distribution from a formula specifying log density.

Construct target distribution from a formula specifying log density.

## Usage

``` r
target_distribution_from_log_density_formula(log_density_formula)
```

## Arguments

- log_density_formula:

  Formula for which right-hand side specifies expression for logarithm
  of (unnormalized) density of target distribution.

## Value

A list with entries

- `log_density`: A function to evaluate log density function for target
  distribution given current position vector.

- `value_and_gradient_log_density`: A function to evaluate value and
  gradient of log density function for target distribution given current
  position vector, returning as a list with entries `value` and
  `gradient`.

## Examples

``` r
target_distribution <- target_distribution_from_log_density_formula(
  ~ (-(x^2 + y^2) / 8 - (x^2 - y)^2 - (x - 1)^2 / 10)
)
target_distribution$value_and_gradient_log_density(c(0.1, -0.3))
#> $value
#> [1] -0.1896
#> 
#> $gradient
#>     x     y 
#> 0.031 0.695 
#> 
```
