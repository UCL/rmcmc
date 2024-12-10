#' Create object to adapt proposal scale to coerce average acceptance rate.
#'
#' @param algorithm String specifying algorithm to use. One of:
#'   * "stochastic_approximation" to use a Robbins-Monro (1951) based scheme,
#'   * "dual_averaging" to use dual-averaging scheme of Nesterov (2009).
#' @param initial_scale Initial value to use for scale parameter. If not set
#'   explicitly a proposal and dimension dependent default will be used.
#' @param target_accept_prob Target value for average accept probability for
#'   chain. If not set a proposal dependent default will be used.
#' @param ... Any additional algorithmic parameters to pass through to
#'   [dual_averaging_scale_adapter()] or [stochastic_approximation_scale_adapter()].
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
#' @seealso [dual_averaging_scale_adapter(), stochastic_approximation_scale_adapter()]
#'
#' @references Nesterov, Y. (2009). Primal-dual subgradient methods for convex
#'   problems. _Mathematical Programming_, 120(1), 221-259.
#' @references Robbins, H., & Monro, S. (1951). A stochastic approximation
#'   method. _The Annals of Mathematical Statistics_, 400-407.
#'
#' @export
#'
#' @examples
#' target_distribution <- list(
#'   log_density = function(x) -sum(x^2) / 2,
#'   grad_log_density = function(x) -x
#' )
#' proposal <- barker_proposal(target_distribution)
#' adapter <- scale_adapter(initial_scale = 1., target_accept_prob = 0.4)
#' adapter$initialize(proposal, chain_state(c(0, 0)))
scale_adapter <- function(
    algorithm = "dual_averaging",
    initial_scale = NULL,
    target_accept_prob = NULL,
    ...) {
  adapter_function <- switch(algorithm,
    dual_averaging = dual_averaging_scale_adapter,
    stochastic_approximation = stochastic_approximation_scale_adapter,
    stop(sprintf("Unrecognized algorithm choice %s"), algorithm)
  )
  adapter_function(initial_scale, target_accept_prob, ...)
}

#' Create object to adapt proposal scale to coerce average acceptance rate using
#' a Robbins and Monro (1951) scheme.
#'
#' When combined with [covariance_shape_adapter()] corresponds to Algorithm 4 in
#' Andrieu and Thoms (2009).
#'
#' @inherit scale_adapter params return
#'
#' @param kappa Decay rate exponent in `[0.5, 1]` for adaptation learning rate.
#'
#' @references Robbins, H., & Monro, S. (1951). A stochastic approximation
#'   method. _The Annals of Mathematical Statistics_, 400-407.
#'
#' @export
#'
#' @examples
#' target_distribution <- list(
#'   log_density = function(x) -sum(x^2) / 2,
#'   grad_log_density = function(x) -x
#' )
#' proposal <- barker_proposal(target_distribution)
#' adapter <- stochastic_approximation_scale_adapter(
#'   initial_scale = 1., target_accept_prob = 0.4
#' )
#' adapter$initialize(proposal, chain_state(c(0, 0)))
stochastic_approximation_scale_adapter <- function(
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
#' @inherit scale_adapter params return
#'
#' @param kappa Decay rate exponent in `[0.5, 1]` for adaptation learning rate.
#'   Defaults to value recommended in Hoffman and Gelman (2014).
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
#' @references Nesterov, Y. (2009). Primal-dual subgradient methods for convex
#'   problems. _Mathematical Programming_, 120(1), 221-259.
#' @references Hoffman, M. D., & Gelman, A. (2014). The No-U-Turn sampler:
#'   adaptively setting path lengths in Hamiltonian Monte Carlo.
#'   _Journal of Machine Learning Research_, 15(1), 1593-1623.
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

#' Create object to adapt proposal shape.
#'
#' @param kappa Decay rate exponent in `[0.5, 1]` for adaptation learning rate.
#'   Value of 1 (default) corresponds to computing empirical variances.
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
#' adapter <- shape_adapter()
#' adapter$initialize(proposal, chain_state(c(0, 0)))
shape_adapter <- function(type = "covariance", kappa = 1) {
  adapter_function <- switch(type,
    covariance = covariance_shape_adapter,
    variance = variance_shape_adapter,
    stop(sprintf("Unrecognized type choice %s"), type)
  )
  adapter_function(kappa)
}


#' Create object to adapt proposal with per dimension scales based on estimates
#' of target distribution variances.
#'
#' Corresponds to variance variant of Algorithm 2 in Andrieu and Thoms (2009),
#' which is itself a restatement of method proposed in Haario et al. (2001).
#'
#' @param kappa Decay rate exponent in `[0.5, 1]` for adaptation learning rate.
#'   Value of 1 (default) corresponds to computing empirical variances.
#'
#' @inherit scale_adapter return
#'
#' @references Andrieu, C., & Thoms, J. (2008). A tutorial on adaptive MCMC.
#'   _Statistics and Computing_, 18, 343-373.
#' @references Haario, H., Saksman, E., & Tamminen, J. (2001). An adaptive
#'   Metropolis algorithm. _Bernoulli_, 7(2): 223-242.
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
variance_shape_adapter <- function(kappa = 1) {
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

#' Create object to adapt proposal with shape based on estimate of target
#' distribution covariance matrix.
#'
#' Corresponds to Algorithm 2 in Andrieu and Thoms (2009), which is itself a
#' restatement of method proposed in Haario et al. (2001).
#'
#' Requires `ramcmc` package to be installed for access to efficient rank-1
#' Cholesky update function [ramcmc::chol_update()].
#'
#' @param kappa Decay rate exponent in `[0.5, 1]` for adaptation learning rate.
#'  Value of 1 (default) corresponds to computing empirical covariance matrix.
#'
#' @inherit scale_adapter return
#'
#' @references Andrieu, C., & Thoms, J. (2008). A tutorial on adaptive MCMC.
#'   _Statistics and Computing_, 18, 343-373.
#' @references Haario, H., Saksman, E., & Tamminen, J. (2001). An adaptive
#'   Metropolis algorithm. _Bernoulli_, 7(2): 223-242.
#'
#' @export
#' @examples
#' target_distribution <- list(
#'   log_density = function(x) -sum(x^2) / 2,
#'   grad_log_density = function(x) -x
#' )
#' proposal <- barker_proposal(target_distribution)
#' adapter <- covariance_shape_adapter()
#' adapter$initialize(proposal, chain_state(c(0, 0)))
covariance_shape_adapter <- function(kappa = 1) {
  rlang::check_installed("ramcmc", reason = "to use this function")
  mean_estimate <- NULL
  chol_covariance_estimate <- NULL
  initialize <- function(proposal, initial_state) {
    mean_estimate <<- initial_state$position()
    chol_covariance_estimate <<- diag(1., initial_state$dimension())
  }
  update <- function(proposal, sample_index, state_and_statistics) {
    # Offset sample_index by 1 so that initial identity covariance estimate acts
    # as regularizer
    beta <- (sample_index + 1)^(-kappa)
    position <- state_and_statistics$state$position()
    mean_estimate <<- mean_estimate + beta * (position - mean_estimate)
    chol_covariance_estimate <<- ramcmc::chol_update(
      sqrt(1 - beta) * chol_covariance_estimate,
      sqrt(beta) * (position - mean_estimate)
    )
    proposal$update(shape = chol_covariance_estimate)
  }
  list(
    initialize = initialize,
    update = update,
    finalize = NULL,
    state = function() {
      list(
        mean_estimate = mean_estimate,
        chol_covariance_estimate = chol_covariance_estimate
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
