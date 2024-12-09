#' Sample new state from Barker proposal.
#'
#' @inheritParams sample_metropolis_hastings
#' @param state Current chain state.
#' @param scale_and_shape Scalar, vector or matrix which scales and shapes
#'  proposal distribution. If a scalar (in which case the value should be
#'  non-negative) the auxiliary vector will be isotropically scaled by the
#'  value. If a vector (in which case the value should be equal in length to the
#'  dimension of the space and all entries non-negative) each dimension of the
#'  auxiliary vector will be scaled separately. If a matrix (in which case the
#'  value should be a square matrix with size equal to the dimension of the
#'  space) then by pre-multiplying the auxiliary vector arbitrary linear
#'  transformations can be performed.
#' @param sample_auxiliary Function which generates a random vector from
#'   auxiliary variable distribution.
#' @param sample_uniform Function which generates a random vector from standard
#'   uniform distribution given an integer size.
#'
#' @return Proposed new chain state.
#'
#' @keywords internal
sample_barker <- function(
    state,
    target_distribution,
    scale_and_shape,
    sample_auxiliary = stats::rnorm,
    sample_uniform = stats::runif) {
  grad <- state$gradient_log_density(target_distribution)
  dim <- state$dimension()
  auxiliary <- sample_auxiliary(dim)
  p_signs <- logistic_sigmoid((grad %@% scale_and_shape) * auxiliary)
  signs <- 2 * (sample_uniform(dim) < p_signs) - 1
  momentum <- signs * auxiliary
  position <- state$position() + (scale_and_shape %@% momentum)
  chain_state(position = position, momentum = momentum)
}

#' Compute logarithm of Barker proposal density ratio.
#'
#' @inheritParams sample_barker
#' @param proposed_state Proposed chain state.
#'
#' @return Logarithm of proposal density ratio.
#'
#' @keywords internal
log_density_ratio_barker <- function(
    state,
    proposed_state,
    target_distribution,
    scale_and_shape) {
  sum(
    log1p_exp(
      state$momentum() * (
        state$gradient_log_density(target_distribution) %@% scale_and_shape
      )
    ) - log1p_exp(
      proposed_state$momentum() * (
        proposed_state$gradient_log_density(target_distribution) %@% scale_and_shape
      )
    )
  )
}

#' Create a new Barker proposal object.
#'
#' Returns a list with function to sample from the proposal, evaluate the log
#' density ratio for a state pair for the proposal and update the proposal
#' parameters. The proposal has two parameters `scale` and `shape`. At least one
#' of `scale` and `shape` must be set before sampling from the proposal or
#' evaluating the log density ratio.
#'
#' @inheritParams sample_barker
#' @param scale Scale parameter of proposal distribution. A non-negative scalar
#'   value determining scale of steps proposed.
#' @param shape Shape parameter of proposal distribution. Either a vector
#'    corresponding to a diagonal shape matrix with per-dimension scaling
#'    factors, or a matrix allowing arbitrary linear transformations.
#'
#' @return Proposal object. A list with entries
#' * `sample`: a function to generate sample from proposal distribution given
#'   current chain state,
#' * `log_density_ratio`: a  function to compute log density ratio for proposal
#'   for a given pair of current and proposed chain states,
#' * `update`: a function to update parameters of proposal,
#' * `parameters`: a function to return list of current parameter values.
#' * `default_target_accept_prob`: a function returning the default target
#'   acceptance rate to use for any scale adaptation.
#' * `default_initial_scale`: a function which given a dimension gives a default
#'   value to use for the initial proposal scale parameter.
#'
#' @export
#'
#' @examples
#' target_distribution <- list(
#'   log_density = function(x) -sum(x^2) / 2,
#'   gradient_log_density = function(x) -x
#' )
#' proposal <- barker_proposal(target_distribution, scale = 1.)
#' state <- chain_state(c(0., 0.))
#' withr::with_seed(876287L, proposed_state <- proposal$sample(state))
#' log_density_ratio <- proposal$log_density_ratio(state, proposed_state)
#' proposal$update(scale = 0.5)
barker_proposal <- function(
    target_distribution,
    scale = NULL,
    shape = NULL,
    sample_auxiliary = stats::rnorm,
    sample_uniform = stats::runif) {
  scale_and_shape_proposal(
    sample = function(state, scale_and_shape) {
      sample_barker(
        state, target_distribution, scale_and_shape, sample_auxiliary, sample_uniform
      )
    },
    log_density_ratio = function(state, proposed_state, scale_and_shape) {
      log_density_ratio_barker(
        state, proposed_state, target_distribution, scale_and_shape
      )
    },
    scale = scale,
    shape = shape,
    default_target_accept_prob = 0.574,
    default_initial_scale = function(dimension) dimension^(-1 / 6)
  )
}
