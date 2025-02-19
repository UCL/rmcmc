% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/barker.R
\name{bimodal_barker_proposal}
\alias{bimodal_barker_proposal}
\title{Create a new Barker proposal object with bimodal noise distribution.}
\usage{
bimodal_barker_proposal(
  sigma = 0.1,
  scale = NULL,
  shape = NULL,
  sample_uniform = stats::runif
)
}
\arguments{
\item{sigma}{Standard deviation of equally-weighted normal components in
bimodal auxiliary noise distribution, with corresponding means of
\verb{±sqrt(1 - sigma^2)}.}

\item{scale}{Scale parameter of proposal distribution. A non-negative scalar
value determining scale of steps proposed.}

\item{shape}{Shape parameter of proposal distribution. Either a vector
corresponding to a diagonal shape matrix with per-dimension scaling
factors, or a matrix allowing arbitrary linear transformations.}

\item{sample_uniform}{Function which generates a random vector from standard
uniform distribution given an integer size.}
}
\value{
Proposal object. A list with entries
\itemize{
\item \code{sample}: a function to generate sample from proposal distribution given
current chain state,
\item \code{log_density_ratio}: a  function to compute log density ratio for proposal
for a given pair of current and proposed chain states,
\item \code{update}: a function to update parameters of proposal,
\item \code{parameters}: a function to return list of current parameter values.
\item \code{default_target_accept_prob}: a function returning the default target
acceptance rate to use for any scale adaptation.
\item \code{default_initial_scale}: a function which given a dimension gives a default
value to use for the initial proposal scale parameter.
}
}
\description{
Convenience function for creating a Barker proposal with bimodal auxiliary
noise variable distribution, corresponding to equally-weighted normal
components with shared variance \code{sigma} and means \verb{±sqrt(1 - sigma^2)}.
This choice of noise distribution was suggested in Vogrinc et al. (2023) and
found to give improved performance over the default choice of a standard
normal auxiliary noise distribution in a range of targets.
}
\details{
For more details see the vignette:
\code{vignette("adjusting-noise-distribution", package = "rmcmc")}
}
\examples{
target_distribution <- list(
  log_density = function(x) -sum(x^2) / 2,
  gradient_log_density = function(x) -x
)
proposal <- bimodal_barker_proposal(scale = 1.)
state <- chain_state(c(0., 0.))
withr::with_seed(
  876287L, proposed_state <- proposal$sample(state, target_distribution)
)
log_density_ratio <- proposal$log_density_ratio(
  state, proposed_state, target_distribution
)
proposal$update(scale = 0.5)
}
\references{
Vogrinc, J., Livingstone, S., & Zanella, G. (2023). Optimal
design of the Barker proposal and other locally balanced
Metropolis–Hastings algorithms. \emph{Biometrika}, 110(3), 579-595.
}
\seealso{
\code{\link[=barker_proposal]{barker_proposal()}}
}
