#' Create object to adapt proposal scale to coerce average acceptance rate using
#' a Robbins and Monro (1951) scheme.
#'
#' @param initial_scale Initial value to use for scale parameter. If not set
#'   explicitly a proposal and dimension dependent default will be used.
#' @param target_accept_prob Target value for average accept probability for
#'   chain. If not set a proposal dependent default will be used.
#' @param kappa Decay rate exponent in `[0.5, 1]` for adaptation learning rate.
#'
#' @return List of functions with entries
#' * `initialize`, a function for initializing adapter state and proposal
#'   parameters at beginning of chain,
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
#' adapter <- simple_scale_adapter(initial_scale = 1., target_accept_prob = 0.4)
#' adapter$initialize(proposal, chain_state(c(0, 0)))
simple_scale_adapter <- function(
    initial_scale = NULL, target_accept_prob = NULL, kappa = 0.6) {
  log_scale <- NULL
  initialize <- function(proposal, initial_state) {
    if (is.null(initial_scale)) {
      initial_scale <- proposal$default_initial_scale(initial_state$dimension())
    }
    log_scale <<- log(initial_scale)
    proposal$update(scale = initial_scale)
  }
  update <- function(proposal, sample_index, state_and_statistics) {
    if (is.null(target_accept_prob)) {
      target_accept_prob <- proposal$default_target_accept_prob()
    }
    beta <- sample_index^(-kappa)
    accept_prob <- state_and_statistics$statistics$accept_prob
    log_scale <<- log_scale + beta * (accept_prob - target_accept_prob)
    proposal$update(scale = exp(log_scale))
  }
  list(
    initialize = initialize,
    update = update,
    finalize = NULL,
    state = function() list(log_scale = log_scale)
  )
}

#' Create object to adapt proposal scale to coerce average acceptance rate
#' using dual averaging scheme of Nesterov (2009) and Hoffman and Gelman (2014).
#'
#' @inherit simple_scale_adapter params return
#'
#' @param gamma Regularization coefficient for (log) scale in dual averaging
#'   algorithm. Controls amount of regularization of (log) scale towards `mu`.
#'   Should be set to a non-negative value. Defaults to value recommended in
#'   Hoffman and Gelman (2014).
#' @param iteration_offset Offset to chain iteration used for the iteration
#'   based weighting of the adaptation statistic error estimate. Should be set
#'   to a non-negative value. A value greater than zero has the effect of
#'   stabilizing early iterations. Defaults to value recommended in
#'   Hoffman and Gelman (2014).
#' @param mu Value to regularize (log) scale towards. If `NULL` (the default),
#'   `mu` will be set to `log(10 * initial_scale)`, as recommended in Hoffman
#'   and Gelman (2014).
#'
#' @export
#'
#' @examples
#' target_distribution <- list(
#'   log_density = function(x) -sum(x^2) / 2,
#'   grad_log_density = function(x) -x
#' )
#' proposal <- barker_proposal(target_distribution)
#' adapter <- dual_averaging_scale_adapter(
#'   initial_scale = 1., target_accept_prob = 0.4
#' )
#' adapter$initialize(proposal, chain_state(c(0, 0)))
dual_averaging_scale_adapter <- function(
    initial_scale = NULL,
    target_accept_prob = NULL,
    kappa = 0.75,
    gamma = 0.05,
    iteration_offset = 10,
    mu = NULL) {
  log_scale <- NULL
  smoothed_log_scale <- 0
  accept_prob_error <- 0
  initialize <- function(proposal, initial_state) {
    if (is.null(initial_scale)) {
      initial_scale <- proposal$default_initial_scale(initial_state$dimension())
    }
    if (is.null(mu)) {
      mu <<- log(10 * initial_scale)
    }
    log_scale <<- log(initial_scale)
    proposal$update(scale = initial_scale)
  }
  update <- function(proposal, sample_index, state_and_statistics) {
    if (is.null(target_accept_prob)) {
      target_accept_prob <- proposal$default_target_accept_prob()
    }
    accept_prob <- state_and_statistics$statistics$accept_prob
    offset_sample_index <- sample_index + iteration_offset
    accept_prob_error <<- (
      (1 - 1 / offset_sample_index) * accept_prob_error + (
        target_accept_prob - accept_prob
      ) / offset_sample_index
    )
    log_scale <<- mu - sqrt(sample_index) * accept_prob_error / gamma
    beta <- sample_index^(-kappa)
    smoothed_log_scale <<- beta * log_scale + (1 - beta) * smoothed_log_scale
    proposal$update(scale = exp(log_scale))
  }
  finalize <- function(proposal) {
    proposal$update(scale = exp(smoothed_log_scale))
  }
  list(
    initialize = initialize,
    update = update,
    finalize = finalize,
    state = function() {
      list(
        log_scale = log_scale,
        smoothed_log_scale = smoothed_log_scale,
        accept_prob_error = accept_prob_error
      )
    }
  )
}

#' Create object to adapt proposal with per dimension scales based on estimates
#' of target distribution variances.
#'
#' @param kappa Decay rate exponent in `[0.5, 1]` for adaptation learning rate.
#'
#' @inherit simple_scale_adapter return
#'
#' @export
#' @examples
#' target_distribution <- list(
#'   log_density = function(x) -sum(x^2) / 2,
#'   grad_log_density = function(x) -x
#' )
#' proposal <- barker_proposal(target_distribution)
#' adapter <- variance_shape_adapter()
#' adapter$initialize(proposal, chain_state(c(0, 0)))
variance_shape_adapter <- function(kappa = 0.6) {
  mean_estimate <- NULL
  variance_estimate <- NULL
  initialize <- function(proposal, initial_state) {
    mean_estimate <<- initial_state$position()
    variance_estimate <<- rep(1., initial_state$dimension())
  }
  update <- function(proposal, sample_index, state_and_statistics) {
    # Offset sample_index by 1 so that initial unity variance_estimate acts as
    # regularizer
    beta <- (sample_index + 1)^(-kappa)
    position <- state_and_statistics$state$position()
    mean_estimate <<- mean_estimate + beta * (position - mean_estimate)
    variance_estimate <<- variance_estimate + beta * (
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
#' @inheritParams simple_scale_adapter
#'
#' @inherit simple_scale_adapter return
#'
#' @export
#'
#' @examples
#' target_distribution <- list(
#'   log_density = function(x) -sum(x^2) / 2,
#'   grad_log_density = function(x) -x
#' )
#' proposal <- barker_proposal(target_distribution)
#' adapter <- robust_shape_adapter(initial_scale = 1., target_accept_prob = 0.4)
#' adapter$initialize(proposal, chain_state(c(0, 0)))
robust_shape_adapter <- function(
    initial_scale = NULL, target_accept_prob = NULL, kappa = 0.6) {
  rlang::check_installed("ramcmc", reason = "to use this function")
  shape <- NULL
  initialize <- function(proposal, initial_state) {
    if (is.null(initial_scale)) {
      initial_scale <- proposal$default_initial_scale(initial_state$dimension())
    }
    shape <<- initial_scale * diag(initial_state$dimension())
    proposal$update(shape = shape)
  }
  update <- function(proposal, sample_index, state_and_statistics) {
    if (is.null(target_accept_prob)) {
      target_accept_prob <- proposal$default_target_accept_prob()
    }
    momentum <- state_and_statistics$proposed_state$momentum()
    accept_prob <- state_and_statistics$statistics$accept_prob
    shape <<- ramcmc::adapt_S(
      shape, momentum, accept_prob, sample_index, target_accept_prob, kappa
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
