#' Construct target distribution from a BridgeStan `StanModel` object.
#'
#' @param model Stan model object to use for target (posterior) distribution.
#'
#' @return A list with entries
#' * `log_density`: A function to evaluate log density function for target
#'   distribution given current position vector.
#' * `value_and_gradient_log_density`: A function to evaluate value and gradient of
#'   log density function for target distribution given current position vector,
#'   returning as a list with entries `value` and `gradient`.
#'
#' @export
#'
#' @examplesIf requireNamespace("bridgestan", quietly = TRUE)
#' model <- example_gaussian_stan_model()
#' target_distribution <- target_distribution_from_stan_model(model)
#' withr::with_seed(
#'   876287L, state <- chain_state(stats::rnorm(model$param_unc_num()))
#' )
#' state$log_density(target_distribution)
target_distribution_from_stan_model <- function(model) {
  list(
    log_density = model$log_density,
    value_and_gradient_log_density = function(position) {
      value_and_gradient <- model$log_density_gradient(position)
      names(value_and_gradient) <- c("value", "gradient")
      value_and_gradient
    }
  )
}

#' Construct trace function from a BridgeStan `StanModel` object.
#'
#' @param model Stan model object to use to generate (constrained) parameters to
#'   trace.
#' @param include_log_density Whether to include an entry `log_density`
#'   corresponding to current log density for target distribution in values
#'   returned by trace function.
#' @param include_generated_quantities Whether to included generated quantities
#'   in Stan model definition in values returned by trace function.
#' @param include_transformed_parameters Whether to include transformed
#'   parameters in Stan model definition in values returned by trace function.
#'
#' @return A function which given `chain_state` object returns a named vector of
#'   values to trace during sampling. The constrained parameter values of model
#'   will always be included.
#'
#' @export
#'
#' @examplesIf requireNamespace("bridgestan", quietly = TRUE)
#' model <- example_gaussian_stan_model()
#' trace_function <- trace_function_from_stan_model(model)
#' withr::with_seed(876287L, state <- chain_state(rnorm(model$param_unc_num())))
#' trace_function(state)
trace_function_from_stan_model <- function(
    model,
    include_log_density = TRUE,
    include_generated_quantities = FALSE,
    include_transformed_parameters = FALSE) {
  function(state) {
    position <- state$position()
    trace_values <- model$param_constrain(
      position, include_transformed_parameters, include_generated_quantities
    )
    names(trace_values) <- model$param_names()
    if (include_log_density) {
      trace_values["log_density"] <- model$log_density(position)
    }
    trace_values
  }
}

#' Construct an example BridgeStan `StanModel` object for a Gaussian model.
#'
#' Requires BridgeStan package to be installed. Generative model is assumed to
#  be of the form `y ~ normal(mu, sigma)` for unknown `mu` and `sigma`.
#'
#' @param n_data Number of independent data points `y` to generate and condition
#'   model against.
#' @param seed Integer seed for Stan model.
#'
#' @return BridgeStan StanModel object.
#'
#' @export
#'
#' @examplesIf requireNamespace("bridgestan", quietly = TRUE)
#' model <- example_gaussian_stan_model(n_data = 5)
#' model$param_names()
example_gaussian_stan_model <- function(n_data = 50, seed = 1234L) {
  rlang::check_installed("bridgestan", reason = "to use this function")
  model_string <- "data {
  int<lower=0> N;
  vector[N] y;
  }
  parameters {
    real mu;
    real<lower=0> sigma;
  }
  model {
    y ~ normal(mu, sigma);
  }"
  withr::with_seed(seed, y <- stats::rnorm(n_data))
  data_string <- sprintf('{"N": %i, "y": [%s]}', n_data, toString(y))
  model_file <- tempfile("gaussian", fileext = ".stan")
  withr::with_tempfile("model_file",
    {
      writeLines(model_string, model_file)
      bridgestan::StanModel$new(model_file, data_string, seed)
    },
    pattern = "gaussian",
    fileext = ".stan"
  )
}
