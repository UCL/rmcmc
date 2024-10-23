#' Create object to adapt proposal scale to coerce average acceptance rate.
#'
#' @param proposal Proposal object to adapt. Must define an `update` function
#'   which accepts a parameter `scale` for setting scale parameter of proposal.
#' @param initial_scale Initial value to use for scale parameter. If not set
#'   explicitly a proposal and dimension dependent default will be used.
#' @param target_accept_prob Target value for average accept probability for
#'   chain. If not set a proposal dependent default will be used.
#' @param kappa Decay rate exponent in `[0.5, 1]` for adaptation learning rate.
#'
#' @return List of functions with entries
#' * `initialize`, a function for initializing adapter state at beginning of
#'   chain,
#' * `update` a function for updating adapter state and proposal parameters on
#'   each chain iteration,
#' * `finalize` a function for performing any final updates to adapter state and
#'   proposal parameters on completion of chain sampling (may be `NULL` if
#'   unused).
#' * `state` a zero-argument function for accessing current values of adapter
#'   state variables.
#'
#' @export
#'
#' @examples
#' target_distribution <- list(
#'   log_density = function(x) -sum(x^2) / 2,
#'   grad_log_density = function(x) -x
#' )
#' proposal <- barker_proposal(target_distribution)
#' adapter <- scale_adapter(
#'   proposal,
#'   initial_scale = 1., target_accept_prob = 0.4
#' )
scale_adapter <- function(
    proposal, initial_scale = NULL, target_accept_prob = NULL, kappa = 0.6) {
  log_scale <- NULL
  if (is.null(target_accept_prob)) {
    target_accept_prob <- proposal$default_target_accept_prob()
  }
  initialize <- function(initial_state) {
    if (is.null(initial_scale)) {
      initial_scale <- proposal$default_initial_scale(initial_state$dimension())
    }
    log_scale <<- log(initial_scale)
    proposal$update(scale = initial_scale)
  }
  update <- function(sample_index, state_and_statistics) {
    gamma <- sample_index^(-kappa)
    accept_prob <- state_and_statistics$statistics$accept_prob
    log_scale <<- log_scale + gamma * (accept_prob - target_accept_prob)
    proposal$update(scale = exp(log_scale))
  }
  list(
    initialize = initialize,
    update = update,
    finalize = function() {},
    state = function() list(log_scale = log_scale)
  )
}

#' Create object to adapt proposal with per dimension scales based on estimates
#' of target distribution variances.
#'
#' @param proposal Proposal object to adapt. Must define an `update` function
#'   which accepts a parameter `shape` for setting shape parameter of proposal.
#' @param kappa Decay rate exponent in `[0.5, 1]` for adaptation learning rate.
#'
#' @inherit scale_adapter return
#'
#' @export
#' @examples
#' target_distribution <- list(
#'   log_density = function(x) -sum(x^2) / 2,
#'   grad_log_density = function(x) -x
#' )
#' proposal <- barker_proposal(target_distribution)
#' adapter <- variance_adapter(proposal)
variance_adapter <- function(proposal, kappa = 0.6) {
  mean_estimate <- NULL
  variance_estimate <- NULL
  initialize <- function(initial_state) {
    mean_estimate <<- initial_state$position()
    variance_estimate <<- rep(1., initial_state$dimension())
  }
  update <- function(sample_index, state_and_statistics) {
    gamma <- sample_index^(-kappa)
    position <- state_and_statistics$state$position()
    mean_estimate <<- mean_estimate + gamma * (position - mean_estimate)
    variance_estimate <<- variance_estimate + gamma * (
      (position - mean_estimate)^2 - variance_estimate
    )
    proposal$update(shape = sqrt(variance_estimate))
  }
  list(
    initialize = initialize,
    update = update,
    finalize = NULL,
    state = function() {
      list(
        mean_estimate = mean_estimate, variance_estimate = variance_estimate
      )
    }
  )
}

#' Create object to adapt proposal shape (and scale) using robust adaptive
#' Metropolis algorithm of Vihola (2012).
#'
#' Requires `ramcmc` package to be installed.
#'
#' @references Vihola, M. (2012). Robust adaptive Metropolis algorithm with
#'     coerced acceptance rate. _Statistics and Computing_, 22, 997-1008.
#'     <https://doi.iorg/10.1007/s11222-011-9269-5>
#'
#' @inheritParams scale_adapter
#'
#' @inherit scale_adapter return
#'
#' @export
#'
#' @examples
#' target_distribution <- list(
#'   log_density = function(x) -sum(x^2) / 2,
#'   grad_log_density = function(x) -x
#' )
#' proposal <- barker_proposal(target_distribution)
#' adapter <- robust_shape_adapter(
#'   proposal,
#'   initial_scale = 1.,
#'   target_accept_prob = 0.4
#' )
robust_shape_adapter <- function(
    proposal, initial_scale, target_accept_prob = 0.4, kappa = 0.6) {
  rlang::check_installed("ramcmc", reason = "to use this function")
  shape <- NULL
  initialize <- function(initial_state) {
    shape <<- initial_scale * diag(initial_state$dimension())
    proposal$update(shape = shape)
  }
  update <- function(sample_index, state_and_statistics) {
    momentum <- state_and_statistics$proposed_state$momentum()
    accept_prob <- state_and_statistics$statistics$accept_prob
    shape <<- ramcmc::adapt_S(
      shape, momentum, accept_prob, sample_index - 1, target_accept_prob, kappa
    )
    proposal$update(shape = shape)
  }
  list(
    initialize = initialize,
    update = update,
    finalize = NULL,
    state = function() list(shape = shape)
  )
}
