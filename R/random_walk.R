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


#' Create a new random walk proposal object.
#'
#' @inherit barker_proposal return params description
#'
#' @export
#'
#' @examples
#' target_distribution <- list(
#'   log_density = function(x) -sum(x^2) / 2
#' )
#' proposal <- random_walk_proposal(target_distribution, scale = 1.)
#' state <- chain_state(c(0., 0.))
#' withr::with_seed(876287L, proposed_state <- proposal$sample(state))
#' log_density_ratio <- proposal$log_density_ratio(state, proposed_state)
#' proposal$update(scale = 0.5)
random_walk_proposal <- function(
    target_distribution,
    scale = NULL,
    shape = NULL,
    sample_auxiliary = stats::rnorm) {
  scale_and_shape_proposal(
    sample = function(state, scale_and_shape) {
      sample_random_walk(
        state, target_distribution, scale_and_shape, sample_auxiliary
      )
    },
    log_density_ratio = function(state, proposed_state, scale_and_shape) 0,
    scale = scale,
    shape = shape,
    default_target_accept_prob = 0.234,
    default_initial_scale = function(dimension) 2.38 / (dimension)^(1 / 2)
  )
}
