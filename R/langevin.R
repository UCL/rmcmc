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
  momentum <- state$momentum() + Matrix::drop(
    Matrix::t(scale_and_shape) %*% gradient
  ) / 2
  position <- state$position() + Matrix::drop(scale_and_shape %*% momentum)
  state <- chain_state(position = position, momentum = momentum)
  gradient <- state$gradient_log_density(target_distribution)
  momentum <- state$momentum() + Matrix::drop(
    Matrix::t(scale_and_shape) %*% gradient
  ) / 2
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
  sum(proposed_state$momentum()^2 - state$momentum()^2) / 2
}

#' Create a new Langevin proposal object.
#'
#' @inherit barker_proposal return params description
#'
#' @export
#'
#' @examples
#' target_distribution <- list(
#'   log_density = function(x) -sum(x^2) / 2,
#'   gradient_log_density = function(x) -x
#' )
#' proposal <- langevin_proposal(target_distribution, scale = 1.)
#' state <- chain_state(c(0., 0.))
#' withr::with_seed(876287L, proposed_state <- proposal$sample(state))
#' log_density_ratio <- proposal$log_density_ratio(state, proposed_state)
#' proposal$update(scale = 0.5)
langevin_proposal <- function(
    target_distribution,
    scale = NULL,
    shape = NULL,
    sample_auxiliary = stats::rnorm) {
  scale_and_shape_proposal(
    sample = function(state, scale_and_shape) {
      sample_langevin(
        state, target_distribution, scale_and_shape, sample_auxiliary
      )
    },
    log_density_ratio = function(state, proposed_state, scale_and_shape) {
      log_density_ratio_langevin(
        state, proposed_state, target_distribution, scale_and_shape
      )
    },
    scale = scale,
    shape = shape
  )
}
