#' Sample new state from random walk proposal.
#'
#' @inherit barker_proposal return params
#'
#' @keywords internal
sample_random_walk <- function(
    state,
    target_distribution,
    scale_and_shape,
    sample_auxiliary = stats::rnorm) {
  momentum <- sample_auxiliary(state$dimension())
  position <- state$position() + (scale_and_shape %@% momentum)
  chain_state(position = position, momentum = momentum)
}

#' Compute logarithm of random_walk proposal density ratio.
#'
#' @inherit log_density_ratio_barker return params
#'
#' @keywords internal
log_density_ratio_random_walk <- function(
    state,
    proposed_state,
    target_distribution,
    scale_and_shape) {
  0
}

#' Create a new (Gaussian) random walk proposal object.
#'
#' The Gaussian random walk proposal samples a new proposed state by perturbing
#' the current state with zero-mean normally distributed noise.
#'
#' @inherit barker_proposal return params
#'
#' @export
#'
#' @examples
#' target_distribution <- list(log_density = function(x) -sum(x^2) / 2)
#' proposal <- random_walk_proposal(scale = 1.)
#' state <- chain_state(c(0., 0.))
#' withr::with_seed(
#'   876287L, proposed_state <- proposal$sample(state, target_distribution)
#' )
#' log_density_ratio <- proposal$log_density_ratio(
#'   state, proposed_state, target_distribution
#' )
#' proposal$update(scale = 0.5)
random_walk_proposal <- function(
    scale = NULL,
    shape = NULL,
    sample_auxiliary = stats::rnorm) {
  scale_and_shape_proposal(
    sample = function(state, target_distribution, scale_and_shape) {
      sample_random_walk(
        state, target_distribution, scale_and_shape, sample_auxiliary
      )
    },
    log_density_ratio = log_density_ratio_random_walk,
    scale = scale,
    shape = shape,
    default_target_accept_prob = 0.234,
    default_initial_scale = function(dimension) 2.38 / (dimension)^(1 / 2)
  )
}
