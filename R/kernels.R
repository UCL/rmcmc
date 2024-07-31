#' Sample from Metropolis-Hastings kernel.
#'
#' @param state Current chain state.
#' @param target_distribution Target stationary distribution for chain. A list
#'  with named entries `log_density` and `gradient_log_density` corresponding to
#'  respectively functions for evaluating the logarithm of the (potentially
#'  unnormalized) density of the target distribution and its gradient.
#'  As an alternative to `gradient_log_density` an entry
#'  `value_and_gradient_log_density` may instead be provided which is a function
#'  returning both the value and gradient of the logarithm of the (unnormalized)
#'  density of the target distribution as a list under the names `value` and
#'  `gradient` respectively.
#' @param proposal Proposal distribution object. Must define entries `sample`, a
#'   function to generate sample from proposal distribution given current chain
#'   state and `log_density_ratio`, a function to compute log density ratio for
#'   proposal for a given pair of current and proposed chain states.
#' @param sample_uniform Function which generates a random vector from standard
#'   uniform distribution given an integer size.
#'
#' @return List with named entries
#' * `state`: corresponding to new chain state,
#' * `proposed_state`: corresponding to proposed chain state,
#' * `statistics`: a list with named entries for statistics of transition, here
#'   this consisting of a named entry `accept_prob` for the Metropolis
#'   acceptance probability.
#'
#' @keywords internal
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
  list(
    state = new_state,
    proposed_state = proposed_state,
    statistics = list(accept_prob = accept_prob))
}
