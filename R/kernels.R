#' Sample from Metropolis-Hastings kernel.
#'
#' @param state Current chain state.
#' @param target_distribution Target distribution for chain.
#' @param proposal Proposal distribution.
#' @param sample_uniform Function which generates a random vector from standard
#'   uniform distribution given an integer size.
#'
#' @return List with named entries `state` corresponding to new chain state and
#'   `statistics`, a list with named entries for statistics of transition, here
#'   this consisting of a named entry `accept_prob` for the Metropolis
#'   acceptance probability.
#' @export
#'
#' @examples
#' target_distribution <- list(
#'   log_density = function(x) sum(x^2) / 2,
#'   grad_log_density = function(x) x
#' )
#' proposal <- barker_proposal(target_distribution, scale = 1.)
#' n_sample <- 1000
#' states <- vector("list", n_sample)
#' states[[1]] <- chain_state(rnorm(2))
#' for (s in 2:n_sample) {
#'   state_and_statistics <- sample_metropolis_hastings(
#'     states[[s - 1]], target_distribution, proposal
#'   )
#'   states[[s]] <- state_and_statistics$state
#' }
sample_metropolis_hastings <- function(
    state,
    target_distribution,
    proposal,
    sample_uniform = stats::runif) {
  proposed_state <- proposal$sample(state)
  log_accept_ratio <- (
    proposed_state$log_density(target_distribution)
    - state$log_density(target_distribution)
    + proposal$log_density_ratio(state, proposed_state)
  )
  accept_prob <- if (is.nan(log_accept_ratio)) 0 else min_1_exp(log_accept_ratio)
  if (sample_uniform(1) < accept_prob) {
    new_state <- proposed_state
  } else {
    new_state <- state
  }
  list(state = new_state, statistics = list(accept_prob = accept_prob))
}
