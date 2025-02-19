% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/kernels.R
\name{sample_metropolis_hastings}
\alias{sample_metropolis_hastings}
\title{Sample from Metropolis-Hastings kernel.}
\usage{
sample_metropolis_hastings(
  state,
  target_distribution,
  proposal,
  sample_uniform = stats::runif
)
}
\arguments{
\item{state}{Current chain state.}

\item{target_distribution}{Target stationary distribution for chain. A list
with named entries \code{log_density} and \code{gradient_log_density} corresponding to
respectively functions for evaluating the logarithm of the (potentially
unnormalized) density of the target distribution and its gradient.
As an alternative to \code{gradient_log_density} an entry
\code{value_and_gradient_log_density} may instead be provided which is a function
returning both the value and gradient of the logarithm of the (unnormalized)
density of the target distribution as a list under the names \code{value} and
\code{gradient} respectively.}

\item{proposal}{Proposal distribution object. Must define entries \code{sample}, a
function to generate sample from proposal distribution given current chain
state and \code{log_density_ratio}, a function to compute log density ratio for
proposal for a given pair of current and proposed chain states.}

\item{sample_uniform}{Function which generates a random vector from standard
uniform distribution given an integer size.}
}
\value{
List with named entries
\itemize{
\item \code{state}: corresponding to new chain state,
\item \code{proposed_state}: corresponding to proposed chain state,
\item \code{statistics}: a list with named entries for statistics of transition, here
this consisting of a named entry \code{accept_prob} for the Metropolis
acceptance probability.
}
}
\description{
Sample from Metropolis-Hastings kernel.
}
\keyword{internal}
