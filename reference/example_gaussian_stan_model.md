# Construct an example BridgeStan `StanModel` object for a Gaussian model.

Requires BridgeStan package to be installed. Generative model is assumed
to be of the form `y ~ normal(mu, sigma)` for unknown
`mu ~ normal(0, 3)` and `sigma ~ half_normal(0, 3)`.

## Usage

``` r
example_gaussian_stan_model(n_data = 50, seed = 1234L)
```

## Arguments

- n_data:

  Number of independent data points `y` to generate and condition model
  against from `normal(0, 1)`.

- seed:

  Integer seed for Stan model.

## Value

BridgeStan StanModel object.

## Examples

``` r
model <- example_gaussian_stan_model(n_data = 5)
#> [1] "BridgeStan not found at location specified by $BRIDGESTAN environment variable, downloading version 2.7.0 to /home/runner/.bridgestan/bridgestan-2.7.0"
#> [1] "Done!"
model$param_names()
#> [1] "mu"    "sigma"
```
