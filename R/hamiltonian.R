#' Sample new state from Hamiltonian proposal.
#'
#' @keywords internal
sample_hamiltonian <- function(
    state,
    target_distribution,
    n_step,
    scale_and_shape,
    sample_auxiliary,
    sample_n_step) {
  state$update(momentum = sample_auxiliary(state))
  if (!is.null(sample_n_step)) {
    n_step <- sample_n_step(n_step)
  }
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
  momentum <- state$momentum() + (gradient %@% scale_and_shape) / 2
  # Inner full 'leapfrog' steps
  for (step in seq_len(n_step - 1)) {
    position <- state$position() + (scale_and_shape %@% momentum)
    state <- chain_state(position = position, momentum = momentum)
    gradient <- state$gradient_log_density(target_distribution)
    momentum <- state$momentum() + (gradient %@% scale_and_shape)
  }
  # Final half-step
  position <- state$position() + (scale_and_shape %@% momentum)
  state <- chain_state(position = position, momentum = momentum)
  gradient <- state$gradient_log_density(target_distribution)
  momentum <- state$momentum() + (gradient %@% scale_and_shape) / 2
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
#' The Hamiltonian proposal augments the target distribution with normally
#' distributed auxiliary momenta variables and simulates the dynamics for a
#' Hamiltonian function corresponding to the negative logarithm of the density
#' of the resulting joint target distribution using a leapfrog integrator, with
#' the proposed new state being the forward integrate state with momenta negated
#' to ensure reversibility.
#'
#' @inherit barker_proposal return params
#'
#' @param n_step Number of leapfrog steps to simulate Hamiltonian dynamics for
#'   in each proposed move, or parameter passed to function specified by
#'   `sample_n_step` argument if not `NULL`.
#' @param sample_auxiliary A function which samples new values for auxiliary
#'   variables (corresponding to a linear transform of momentum) given current
#'   chain state, leaving their standard normal target distribution invariant.
#'   Defaults to a function sampling independent standard normal random variates
#'   but can be used to implement alternative updates such as partial momentum
#'   refreshment. Function should accept a single argument which is passed the
#'   current chain state.
#' @param sample_n_step Optionally a function which randomly samples number of
#'   leapfrog steps to simulate in each proposed move from some integer-valued
#'   distribution, or `NULL` (the default) to use a fixed deterministic number
#'   of steps as specified by `n_step` argument. If a function it should accept
#'   a single argument which will be passed the value of `n_step` which can
#'   be used to specify parameter(s) of distribution to sample number of steps
#'   from.
#'
#' @references Duane, S., Kennedy, A. D., Pendleton, B. J., & Roweth, D. (1987).
#'   Hybrid Monte Carlo. _Physics Letters B_, 195(2), 216-222.
#' @references Neal, R. M. (2011). MCMC Using Hamiltonian Dynamics. In _Handbook
#'   of Markov Chain Monte Carlo_ (pp. 113-162). Chapman and Hall/CRC.
#'
#' @export
#'
#' @examples
#' target_distribution <- list(
#'   log_density = function(x) -sum(x^2) / 2,
#'   gradient_log_density = function(x) -x
#' )
#'
#' # Proposal with fixed number of leapfrog steps
#' proposal <- hamiltonian_proposal(scale = 1., n_step = 5)
#' state <- chain_state(c(0., 0.))
#' withr::with_seed(
#'   876287L, proposed_state <- proposal$sample(state, target_distribution)
#' )
#' log_density_ratio <- proposal$log_density_ratio(
#'   state, proposed_state, target_distribution
#' )
#' proposal$update(scale = 0.5)
#'
#' # Proposal with number of steps randomly sampled uniformly from 5:10
#' sample_uniform_int <- function(lower, upper) {
#'   lower + sample.int(upper - lower + 1, 1) - 1
#' }
#' proposal <- hamiltonian_proposal(
#'   scale = 1.,
#'   n_step = c(5, 10),
#'   sample_n_step = function(n_step) sample_uniform_int(n_step[1], n_step[2])
#' )
#' withr::with_seed(
#'   876287L, proposed_state <- proposal$sample(state, target_distribution)
#' )
#'
#' # Proposal with partial momentum refreshment
#' partial_momentum_update <- function(state, phi = pi / 4) {
#'   momentum <- state$momentum()
#'   if (is.null(momentum)) {
#'     stats::rnorm(state$dimension())
#'   } else {
#'     cos(phi) * momentum + sin(phi) * stats::rnorm(length(momentum))
#'   }
#' }
#' proposal <- hamiltonian_proposal(
#'   scale = 1.,
#'   n_step = 1,
#'   sample_auxiliary = partial_momentum_update
#' )
#' withr::with_seed(
#'   876287L, proposed_state <- proposal$sample(state, target_distribution)
#' )
hamiltonian_proposal <- function(
    n_step,
    scale = NULL,
    shape = NULL,
    sample_auxiliary = function(state) stats::rnorm(state$dimension()),
    sample_n_step = NULL) {
  scale_and_shape_proposal(
    sample = function(state, target_distribution, scale_and_shape) {
      sample_hamiltonian(
        state,
        target_distribution,
        n_step,
        scale_and_shape,
        sample_auxiliary,
        sample_n_step
      )
    },
    log_density_ratio = log_density_ratio_hamiltonian,
    scale = scale,
    shape = shape,
    default_target_accept_prob = 0.8,
    default_initial_scale = function(dimension) dimension^(-1 / 8)
  )
}
