# rmcmc (development version)

# rmcmc 0.1.2

* In `scale_adapter()` and `shape_adapter()`, setting `algorithm` and `type`
  arguments to invalid values now gives informative error message (#88).
* In `target_distribution_from_stan_model()`, `include_generated_quantities`
  argument now works as expected when set to `TRUE` (#94).
* Documentation website now includes article on interfacing with Stan models
  via BridgeStan (#90).

# rmcmc 0.1.1

* Resubmission to CRAN to fix minor issues with formatting of DESCRIPTION and
  not restoring graphical parameters in function in vignette.

# rmcmc 0.1.0

* Initial CRAN submission.
