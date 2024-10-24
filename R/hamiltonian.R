#' Sample new state from Hamiltonian proposal.
#'
#' @keywords internal
sample_hamiltonian <- function(
    state,
    target_distribution,
    n_step,
    scale_and_shape,
    sample_auxiliary = stats::rnorm) {
  state$update(momentum = sample_auxiliary(state$dimension()))
  involution_hamiltonian(state, n_step, scale_and_shape, target_distribution)
}

#' Apply involution underlying Hamiltonian proposal to a chain state.
#'
#' @inheritParams sample_hamiltonian
#'
#' @return Chain state after involution is applied.
#'
#' @keywords internal
involution_hamiltonian <- function(
    state, n_step, scale_and_shape, target_distribution) {
  # Initial half-step
  gradient <- state$gradient_log_density(target_distribution)
  momentum <- state$momentum() + Matrix::drop(
    Matrix::t(scale_and_shape) %*% gradient
  ) / 2
  # Inner full 'leapfrog' steps
  for (step in seq_len(n_step - 1)) {
    position <- state$position() + Matrix::drop(scale_and_shape %*% momentum)
    state <- chain_state(position = position, momentum = momentum)
    gradient <- state$gradient_log_density(target_distribution)
    momentum <- state$momentum() + Matrix::drop(
      Matrix::t(scale_and_shape) %*% gradient
    )
  }
  # Final half-step
  position <- state$position() + Matrix::drop(scale_and_shape %*% momentum)
  state <- chain_state(position = position, momentum = momentum)
  gradient <- state$gradient_log_density(target_distribution)
  momentum <- state$momentum() + Matrix::drop(
    Matrix::t(scale_and_shape) %*% gradient
  ) / 2
  # Negate momentum to make an involution
  state$update(momentum = -momentum)
  state
}

#' Compute logarithm of Hamiltonian proposal density ratio.
#'
#' @inherit log_density_ratio_barker return params
#'
#' @keywords internal
log_density_ratio_hamiltonian <- function(
    state,
    proposed_state,
    target_distribution,
    scale_and_shape) {
  sum(state$momentum()^2 - proposed_state$momentum()^2) / 2
}

#' Create a new Hamiltonian proposal object.
#'
#' @inherit barker_proposal return params description
#'
#' @param n_step Number of leapfrog steps to simulate Hamiltonian dynamics for
#'   in each proposed move.
#'
#' @export
#'
#' @examples
#' target_distribution <- list(
#'   log_density = function(x) -sum(x^2) / 2,
#'   gradient_log_density = function(x) -x
#' )
#' proposal <- hamiltonian_proposal(target_distribution, scale = 1., n_step = 5)
#' state <- chain_state(c(0., 0.))
#' withr::with_seed(876287L, proposed_state <- proposal$sample(state))
#' log_density_ratio <- proposal$log_density_ratio(state, proposed_state)
#' proposal$update(scale = 0.5)
hamiltonian_proposal <- function(
    target_distribution,
    n_step,
    scale = NULL,
    shape = NULL,
    sample_auxiliary = stats::rnorm) {
  scale_and_shape_proposal(
    sample = function(state, scale_and_shape) {
      sample_hamiltonian(
        state, target_distribution, n_step, scale_and_shape, sample_auxiliary
      )
    },
    log_density_ratio = function(state, proposed_state, scale_and_shape) {
      log_density_ratio_hamiltonian(
        state, proposed_state, target_distribution, scale_and_shape
      )
    },
    scale = scale,
    shape = shape,
    default_target_accept_prob = 0.8,
    default_initial_scale = function(dimension) 2.38 / (dimension)^(1 / 4)
  )
}
