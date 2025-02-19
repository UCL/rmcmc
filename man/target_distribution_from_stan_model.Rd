% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bridges.R
\name{target_distribution_from_stan_model}
\alias{target_distribution_from_stan_model}
\title{Construct target distribution from a BridgeStan \code{StanModel} object.}
\usage{
target_distribution_from_stan_model(
  model,
  include_log_density = TRUE,
  include_generated_quantities = FALSE,
  include_transformed_parameters = FALSE
)
}
\arguments{
\item{model}{Stan model object to use for target (posterior) distribution.}

\item{include_log_density}{Whether to include an entry \code{log_density}
corresponding to current log density for target distribution in values
returned by trace function.}

\item{include_generated_quantities}{Whether to included generated quantities
in Stan model definition in values returned by trace function.}

\item{include_transformed_parameters}{Whether to include transformed
parameters in Stan model definition in values returned by trace function.}
}
\value{
A list with entries
\itemize{
\item \code{log_density}: A function to evaluate log density function for target
distribution given current position vector.
\item \code{value_and_gradient_log_density}: A function to evaluate value and gradient
of log density function for target distribution given current position
vector, returning as a list with entries \code{value} and \code{gradient}.
\item \code{trace_function}: A function which given a \code{chain_state()} object returns a
named vector of values to trace during sampling. The constrained parameter
values of model will always be included.
}
}
\description{
Construct target distribution from a BridgeStan \code{StanModel} object.
}
\examples{
\dontshow{if (requireNamespace("bridgestan", quietly = TRUE)) (if (getRversion() >= "3.4") withAutoprint else force)(\{ # examplesIf}
model <- example_gaussian_stan_model()
target_distribution <- target_distribution_from_stan_model(model)
withr::with_seed(
  876287L, state <- chain_state(stats::rnorm(model$param_unc_num()))
)
state$log_density(target_distribution)
target_distribution$trace_function(state)
\dontshow{\}) # examplesIf}
}
