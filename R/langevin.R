#' Sample new state from Langevin proposal.
#'
#' @inherit barker_proposal return params
#'
#' @keywords internal
sample_langevin <- function(
    state,
    target_distribution,
    scale_and_shape,
    sample_auxiliary = stats::rnorm) {
  state$update(momentum = sample_auxiliary(state$dimension()))
  involution_langevin(state, scale_and_shape, target_distribution)
}

#' Apply involution underlying Langevin proposal to a chain state.
#'
#' @inheritParams sample_langevin
#'
#' @return Chain state after involution is applied.
#'
#' @keywords internal
involution_langevin <- function(state, scale_and_shape, target_distribution) {
  gradient <- state$gradient_log_density(target_distribution)
  momentum <- state$momentum() + (gradient %@% scale_and_shape) / 2
  position <- state$position() + (scale_and_shape %@% momentum)
  state <- chain_state(position = position, momentum = momentum)
  gradient <- state$gradient_log_density(target_distribution)
  momentum <- state$momentum() + (gradient %@% scale_and_shape) / 2
  state$update(momentum = -momentum)
  state
}

#' Compute logarithm of Langevin proposal density ratio.
#'
#' @inherit log_density_ratio_barker return params
#'
#' @keywords internal
log_density_ratio_langevin <- function(
    state,
    proposed_state,
    target_distribution,
    scale_and_shape) {
  sum(state$momentum()^2 - proposed_state$momentum()^2) / 2
}

#' Create a new Langevin proposal object.
#'
#' The Langevin proposal is a gradient-based proposal corresponding to a
#' Euler-Maruyama time discretisation of a Langevin diffusion.
#'
#' @inherit barker_proposal return params
#'
#' @references Besag, J. (1994). "Comments on "Representations of knowledge in
#'   complex systems" by U. Grenander and MI Miller". _Journal of the Royal
#'   Statistical Society, Series B_. 56: 591–592.
#' @references Roberts, G. O., & Tweedie, R. L. (1996). Exponential convergence
#'   of Langevin distributions and their discrete approximations. _Bernoulli_ 2
#'   (4), 341 - 363.
#'
#' @export
#'
#' @examples
#' target_distribution <- list(
#'   log_density = function(x) -sum(x^2) / 2,
#'   gradient_log_density = function(x) -x
#' )
#' proposal <- langevin_proposal(scale = 1.)
#' state <- chain_state(c(0., 0.))
#' withr::with_seed(
#'   876287L, proposed_state <- proposal$sample(state, target_distribution)
#' )
#' log_density_ratio <- proposal$log_density_ratio(
#'   state, proposed_state, target_distribution
#' )
#' proposal$update(scale = 0.5)
langevin_proposal <- function(
    scale = NULL,
    shape = NULL,
    sample_auxiliary = stats::rnorm) {
  scale_and_shape_proposal(
    sample = function(state, target_distribution, scale_and_shape) {
      sample_langevin(
        state, target_distribution, scale_and_shape, sample_auxiliary
      )
    },
    log_density_ratio = log_density_ratio_langevin,
    scale = scale,
    shape = shape,
    default_target_accept_prob = 0.574,
    default_initial_scale = function(dimension) dimension^(-1 / 6)
  )
}
