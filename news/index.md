# Changelog

## rmcmc (development version)

## rmcmc 0.1.2

CRAN release: 2025-10-17

- In
  [`scale_adapter()`](http://github-pages.ucl.ac.uk/rmcmc/reference/scale_adapter.md)
  and
  [`shape_adapter()`](http://github-pages.ucl.ac.uk/rmcmc/reference/shape_adapter.md),
  setting `algorithm` and `type` arguments to invalid values now gives
  informative error message
  ([\#88](https://github.com/UCL/rmcmc/issues/88)).
- In
  [`target_distribution_from_stan_model()`](http://github-pages.ucl.ac.uk/rmcmc/reference/target_distribution_from_stan_model.md),
  `include_generated_quantities` argument now works as expected when set
  to `TRUE` ([\#94](https://github.com/UCL/rmcmc/issues/94)).
- Documentation website now includes article on interfacing with Stan
  models via BridgeStan
  ([\#90](https://github.com/UCL/rmcmc/issues/90)).

## rmcmc 0.1.1

CRAN release: 2025-02-04

- Resubmission to CRAN to fix minor issues with formatting of
  DESCRIPTION and not restoring graphical parameters in function in
  vignette.

## rmcmc 0.1.0

- Initial CRAN submission.
