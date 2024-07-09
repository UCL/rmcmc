#' Sample new state from Barker proposal.
#'
#' @param state Current chain state.
#' @param target_distribution Target stationary distribution for chain. A list
#'  with named entries `log_density` and `grad_log_density`corresponding to
#'  respectively functions for evaluating the logarithm of the (potentially
#'  unnormalized) density of the target distribution and its gradient.
#'  As an alternative to `grad_log_density` an entry
#'  `value_and_grad_log_density` may instead be provided which is a function
#'  returning both the value and gradient of the logarithm of the (unnormalized)
#'  density of the target distribution as a list under the names `value` and
#'  `grad` respectively.
#' @param step_size Step size parameter of proposal distribution. A scalar value
#'   determining scale of steps proposed.
#' @param sample_auxiliary Function which generates a random vector from
#'   auxiliary variable distribution.
#' @param sample_uniform Function which generates a random vector from standard
#'   uniform distribution given an integer size.
#'
#' @return Proposed new chain state.
#' @export
#'
#' @examples
#' state <- chain_state(c(0., 0.))
#' target_distribution <- list(
#'   log_density = function(x) sum(x^2) / 2,
#'   grad_log_density = function(x) x
#' )
#' sample_barker(state, target_distribution, step_size = 1.)
sample_barker <- function(
    state,
    target_distribution,
    step_size,
    sample_auxiliary = stats::rnorm,
    sample_uniform = stats::runif) {
  grad <- state$grad_log_density(target_distribution)
  dim <- state$dimension()
  auxiliary <- step_size * sample_auxiliary(dim)
  signs <- 2 * (sample_uniform(dim) < logistic_sigmoid(grad * auxiliary)) - 1
  chain_state(state$position() + signs * auxiliary)
}

#' Compute logarithm of Barker proposal density ratio.
#'
#' @param state Current chain state.
#' @param proposed_state Proposed chain state.
#' @param target_distribution Target stationary distribution for chain.
#'
#' @return Logarithm of proposal density ratio.
#' @export
#'
#' @examples
#' state <- chain_state(c(0., 0.))
#' target_distribution <- list(
#'   log_density = function(x) sum(x^2) / 2,
#'   grad_log_density = function(x) x
#' )
#' proposed_state <- sample_barker(state, target_distribution, step_size = 1.)
#' log_density_ratio_barker(
#'   state, proposed_state, target_distribution
#' )
log_density_ratio_barker <- function(
    state,
    proposed_state,
    target_distribution) {
  sum(
    log1p_exp(
      (state$position() - proposed_state$position())
      * state$grad_log_density(target_distribution)
    ) - log1p_exp(
      (proposed_state$position() - state$position())
      * proposed_state$grad_log_density(target_distribution)
    )
  )
}

#' Create a new Barker proposal object.
#'
#' @param target_distribution Target stationary distribution for chain. A list
#'  with named entries `log_density` and `grad_log_density`corresponding to
#'  respectively functions for evaluating the logarithm of the (potentially
#'  unnormalized) density of the target distribution and its gradient.
#'  As an alternative to `grad_log_density` an entry
#'  `value_and_grad_log_density` may instead be provided which is a function
#'  returning both the value and gradient of the logarithm of the (unnormalized)
#'  density of the target distribution as a list under the names `value` and
#'  `grad` respectively.
#' @param step_size Step size parameter of proposal distribution. A scalar value
#'   determining scale of steps proposed.
#' @param sample_auxiliary Function which generates a random vector from
#'   auxiliary variable distribution.
#' @param sample_uniform Function which generates a random vector from standard
#'   uniform distribution given an integer size.
#'
#' @return List with entries `sample`, a function to generate sample from
#'   proposal distribution given current chain state, `log_density_ratio`, a
#'   function to compute log density ratio for proposal for a given pair of
#'   current and proposed chain states and `update`, a function to update
#'   step size parameter of proposal.
#' @export
#'
#' @examples
#' target_distribution <- list(
#'   log_density = function(x) sum(x^2) / 2,
#'   grad_log_density = function(x) x
#' )
#' proposal <- barker_proposal(target_distribution, step_size = 1.)
#' state <- chain_state(c(0., 0.))
#' proposed_state <- proposal$sample(state)
#' log_density_ratio <- proposal$log_density_ratio(state, proposed_state)
#' proposal$update(step_size = 0.5)
barker_proposal <- function(
    target_distribution,
    step_size,
    sample_auxiliary = stats::rnorm,
    sample_uniform = stats::runif) {
  list(
    sample = function(state) sample_barker(
      state, target_distribution, step_size, sample_auxiliary, sample_uniform
    ),
    log_density_ratio = function(state, proposed_state) log_density_ratio_barker(
      state, proposed_state, target_distribution
    ),
    update = function(step_size) step_size <<- step_size
  )
}
