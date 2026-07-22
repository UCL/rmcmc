#' Create object to adapt proposal scale to coerce average acceptance rate.
#'
#' @param algorithm String specifying algorithm to use. One of:
#'   * "stochastic_approximation" to use a Robbins-Monro (1951) based scheme,
#'   * "dual_averaging" to use dual-averaging scheme of Nesterov (2009).
#'   * "adam" to use the Adam optimizer of Kingma and Ba (2014) applied to the
#'     acceptance-rate residual, following the implementation in the `walnuts`
#'     library.
#' @param initial_scale Initial value to use for scale parameter. If not set
#'   explicitly a proposal and dimension dependent default will be used.
#' @param target_accept_prob Target value for average accept probability for
#'   chain. If not set a proposal dependent default will be used.
#' @param ... Any additional algorithmic parameters to pass through to the
#'   selected adapter constructor: see [dual_averaging_scale_adapter()],
#'   [stochastic_approximation_scale_adapter()] or [adam_scale_adapter()] for
#'   the full list of parameters accepted by each. In practice, most users
#'   tuning the Adam adapter only need to adjust `learning_rate`; the
#'   moment-decay parameters `beta_1`, `beta_2`, `epsilon` and
#'   `learn_rate_decay` have sensible defaults that rarely need adjustment.
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
#' @seealso [dual_averaging_scale_adapter()],
#'   [stochastic_approximation_scale_adapter()], [adam_scale_adapter()]
#'
#' @references Nesterov, Y. (2009). Primal-dual subgradient methods for convex
#'   problems. _Mathematical Programming_, 120(1), 221-259.
#' @references Robbins, H., & Monro, S. (1951). A stochastic approximation
#'   method. _The Annals of Mathematical Statistics_, 400-407.
#' @references Kingma, D. P., & Ba, J. (2014). Adam: A method for stochastic
#'   optimization. _arXiv preprint_ arXiv:1412.6980.
#'
#' @export
#'
#' @examples
#' proposal <- barker_proposal()
#' adapter <- scale_adapter(initial_scale = 1., target_accept_prob = 0.4)
#' adapter$initialize(proposal, chain_state(c(0, 0)))
scale_adapter <- function(
  algorithm = "dual_averaging",
  initial_scale = NULL,
  target_accept_prob = NULL,
  ...
) {
  adapter_function <- switch(algorithm,
    dual_averaging = dual_averaging_scale_adapter,
    stochastic_approximation = stochastic_approximation_scale_adapter,
    adam = adam_scale_adapter,
    stop(sprintf("Unrecognized algorithm choice %s", algorithm))
  )
  adapter_function(initial_scale, target_accept_prob, ...)
}

get_initial_scale <- function(proposal, initial_state) {
  # Prefer the current proposal scale (e.g. carried over from a previous
  # warm-up stage) over the proposal/dimension-dependent default.
  current_scale <- proposal$parameters()$scale
  if (!is.null(current_scale)) {
    current_scale
  } else {
    proposal$default_initial_scale(initial_state$dimension())
  }
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
#' @references Andrieu, C., & Thoms, J. (2008). A tutorial on adaptive MCMC.
#'   _Statistics and Computing_, 18, 343-373.
#' @references Robbins, H., & Monro, S. (1951). A stochastic approximation
#'   method. _The Annals of Mathematical Statistics_, 400-407.
#'
#' @export
#'
#' @examples
#' proposal <- barker_proposal()
#' adapter <- stochastic_approximation_scale_adapter(
#'   initial_scale = 1., target_accept_prob = 0.4
#' )
#' adapter$initialize(proposal, chain_state(c(0, 0)))
stochastic_approximation_scale_adapter <- function(
  initial_scale = NULL, target_accept_prob = NULL, kappa = 0.6
) {
  log_scale <- NULL
  initialize <- function(proposal, initial_state) {
    if (is.null(initial_scale)) {
      initial_scale <- get_initial_scale(proposal, initial_state)
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
#' proposal <- barker_proposal()
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
  mu = NULL
) {
  log_scale <- NULL
  smoothed_log_scale <- 0
  accept_prob_error <- 0
  initialize <- function(proposal, initial_state) {
    if (is.null(initial_scale)) {
      initial_scale <- get_initial_scale(proposal, initial_state)
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

#' Create object to adapt proposal scale to coerce average acceptance rate
#' using the Adam optimizer of Kingma and Ba (2014).
#'
#' Applies the Adam stochastic-gradient optimizer to the log of the proposal
#' scale, using the acceptance-rate residual `target_accept_prob - accept_prob`
#' as the (stochastic) gradient. This corresponds to treating
#' `-0.5 * (accept_prob - target_accept_prob)^2` as the objective and is the
#' same gradient signal used by [dual_averaging_scale_adapter()] and
#' [stochastic_approximation_scale_adapter()], but plugged into the Adam
#' update rule instead of a dual-averaging or Robbins-Monro schedule. Follows
#' the reference implementation in the `walnuts` library (see references
#' below).
#'
#' To match the stability provided by dual averaging, a learning-rate decay
#' `learning_rate / t^learn_rate_decay` is applied to the per-iteration Adam
#' step; a value of `learn_rate_decay = 0` recovers standard Adam, while
#' `learn_rate_decay = 0.5` (the default here) is recommended by the `walnuts`
#' authors for stable convergence.
#'
#' @inherit scale_adapter params return
#'
#' @param learning_rate Learning rate for Adam optimizer (the `alpha`
#'   hyperparameter in the original Adam paper). Controls the magnitude of the
#'   update applied to the log-scale on each iteration. Should be positive.
#'   Defaults to `0.05`, which is a practical setting that gives fast
#'   convergence of the acceptance-rate coercion within typical MCMC warm-up
#'   lengths; the original Adam paper uses `1e-3`.
#' @param beta_1 Exponential decay rate for the first-moment estimate of the
#'   gradient (the `beta_1` hyperparameter in the original Adam paper). Should
#'   be in `[0, 1)`. Defaults to `0.9`.
#' @param beta_2 Exponential decay rate for the second-moment (squared
#'   gradient) estimate (the `beta_2` hyperparameter in the original Adam
#'   paper). Should be in `[0, 1)`. Defaults to `0.999`.
#' @param epsilon Small positive constant added to the square root of the
#'   second-moment estimate in the denominator of the Adam update for
#'   numerical stability. Defaults to `1e-8`.
#' @param learn_rate_decay Exponent controlling the decay of the effective
#'   learning rate across iterations: on iteration `t` the effective learning
#'   rate is `learning_rate / t^learn_rate_decay`. Should be in `[0, 1]`. A
#'   value of `0` recovers standard Adam; the default of `0.5` matches the
#'   value recommended in the `walnuts` implementation.
#'
#' @references Kingma, D. P., & Ba, J. (2014). Adam: A method for stochastic
#'   optimization. _arXiv preprint_ arXiv:1412.6980.
#' @references Reference implementation in `walnuts`:
#'   <https://github.com/flatironinstitute/walnuts/blob/main/include/walnuts/adam.hpp>
#'
#' @export
#'
#' @examples
#' proposal <- barker_proposal()
#' adapter <- adam_scale_adapter(
#'   initial_scale = 1., target_accept_prob = 0.4
#' )
#' adapter$initialize(proposal, chain_state(c(0, 0)))
adam_scale_adapter <- function(
  initial_scale = NULL,
  target_accept_prob = NULL,
  learning_rate = 0.05,
  beta_1 = 0.9,
  beta_2 = 0.999,
  epsilon = 1e-8,
  learn_rate_decay = 0.5
) {
  log_scale <- NULL
  m <- 0
  v <- 0
  beta_1_pow <- 1
  beta_2_pow <- 1
  initialize <- function(proposal, initial_state) {
    if (is.null(initial_scale)) {
      initial_scale <- proposal$default_initial_scale(initial_state$dimension())
    }
    log_scale <<- log(initial_scale)
    m <<- 0
    v <<- 0
    beta_1_pow <<- 1
    beta_2_pow <<- 1
    proposal$update(scale = initial_scale)
  }
  update <- function(proposal, sample_index, state_and_statistics) {
    if (is.null(target_accept_prob)) {
      target_accept_prob <- proposal$default_target_accept_prob()
    }
    accept_prob <- state_and_statistics$statistics$accept_prob
    # Gradient of the squared-error objective in log-scale space. Note the
    # sign: `grad = target - observed` combined with the `log_scale -= ...`
    # update below reproduces the "observed > target => scale increases"
    # behaviour of the other scale adapters.
    grad <- target_accept_prob - accept_prob
    beta_1_pow <<- beta_1_pow * beta_1
    beta_2_pow <<- beta_2_pow * beta_2
    m <<- beta_1 * m + (1 - beta_1) * grad
    v <<- beta_2 * v + (1 - beta_2) * grad^2
    m_hat <- m / (1 - beta_1_pow)
    v_hat <- v / (1 - beta_2_pow)
    effective_learning_rate <- learning_rate / sample_index^learn_rate_decay
    log_scale <<- log_scale - (
      effective_learning_rate * m_hat / (sqrt(v_hat) + epsilon)
    )
    proposal$update(scale = exp(log_scale))
  }
  list(
    initialize = initialize,
    update = update,
    finalize = NULL,
    state = function() list(log_scale = log_scale, m = m, v = v)
  )
}

#' Create object to adapt proposal shape.
#'
#' @param type Type of shape adapter to use. One of:
#'   * "variance": Diagonal shape matrix adaptation based on estimates of target
#'     distribution variances (see [variance_shape_adapter()]),
#'   * "covariance": Dense shape matrix adaptation based on estimates of target
#'     distribution covariance matrix (see [covariance_shape_adapter()]).
#' @param kappa Decay rate exponent in `[0.5, 1]` for adaptation learning rate.
#'   Value of 1 (default) corresponds to computing empirical (co)variances.
#'
#' @seealso [variance_shape_adapter()], [covariance_shape_adapter()]
#'
#' @inherit scale_adapter return
#'
#' @export
#' @examples
#' proposal <- barker_proposal()
#' adapter <- shape_adapter()
#' adapter$initialize(proposal, chain_state(c(0, 0)))
shape_adapter <- function(type = "covariance", kappa = 1) {
  adapter_function <- switch(type,
    covariance = covariance_shape_adapter,
    variance = variance_shape_adapter,
    stop(sprintf("Unrecognized type choice %s", type))
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
#' @param initial_shape Optional numeric vector of length equal to the target
#'   distribution dimension, specifying the per-dimension proposal scales to
#'   use as the initial variance estimate. When supplied, takes precedence over
#'   both any current proposal shape and the default unit initialisation. When
#'   `NULL` (default), the adapter reads the current proposal shape at
#'   initialisation time (to carry over state from a previous warm-up stage)
#'   and falls back to unit variances if no current shape is available.
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
#' proposal <- barker_proposal()
#' adapter <- variance_shape_adapter()
#' adapter$initialize(proposal, chain_state(c(0, 0)))
variance_shape_adapter <- function(kappa = 1, initial_shape = NULL) {
  mean_estimate <- NULL
  variance_estimate <- NULL
  initialize <- function(proposal, initial_state) {
    mean_estimate <<- initial_state$position()
    dim <- initial_state$dimension()
    if (!is.null(initial_shape)) {
      # Priority 1: explicit user-supplied starting shape.
      variance_estimate <<- initial_shape^2
    } else {
      # Priority 2: current proposal shape, if it is a vector of the right
      # length (i.e. from a previous variance_shape_adapter stage).
      # The length(current_shape) == dim guard is important:
      # if the previous stage used a covariance_shape_adapter, the
      # proposal's shape is a matrix, not a vector.
      # Squaring a matrix is meaningless here so we fall back to identity.
      # Priority 3: fall back to unit variances.
      current_shape <- proposal$parameters()$shape
      current_shape_is_vector <- length(current_shape) == dim
      variance_estimate <<- if (!is.null(current_shape) && current_shape_is_vector) {
        current_shape^2
      } else {
        rep(1., dim)
      }
    }
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
#' Cholesky update function `ramcmc::chol_update`.
#'
#' @param kappa Decay rate exponent in `[0.5, 1]` for adaptation learning rate.
#'  Value of 1 (default) corresponds to computing empirical covariance matrix.
#' @param initial_shape Optional lower-triangular matrix with the same
#'   dimensions as the target distribution, specifying the Cholesky factor of
#'   the proposal covariance to use as the initial estimate. When supplied,
#'   takes precedence over both any current proposal shape and the default
#'   identity initialisation. When `NULL` (default), the adapter reads the
#'   current proposal shape at initialisation time (to carry over state from a
#'   previous warm-up stage) and falls back to the identity matrix if no
#'   current shape is available.
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
#' proposal <- barker_proposal()
#' adapter <- covariance_shape_adapter()
#' adapter$initialize(proposal, chain_state(c(0, 0)))
covariance_shape_adapter <- function(kappa = 1, initial_shape = NULL) {
  rlang::check_installed("ramcmc", reason = "to use this function")
  mean_estimate <- NULL
  chol_covariance_estimate <- NULL
  initialize <- function(proposal, initial_state) {
    mean_estimate <<- initial_state$position()
    dim <- initial_state$dimension()
    if (!is.null(initial_shape)) {
      # Priority 1: explicit user-supplied starting Cholesky factor.
      chol_covariance_estimate <<- initial_shape
    } else {
      # Priority 2: current proposal shape, if it is a square matrix of the
      # right dimension (i.e. from a previous shape adapter stage).
      # Priority 3: fall back to identity matrix.
      current_shape <- proposal$parameters()$shape
      current_shape_valid <- is.matrix(current_shape) && nrow(current_shape) == dim
      chol_covariance_estimate <<- if (!is.null(current_shape) && current_shape_valid) {
        current_shape
      } else {
        diag(1., dim)
      }
    }
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
#'     \doi{10.1007/s11222-011-9269-5}
#'
#' @inheritParams stochastic_approximation_scale_adapter
#'
#' @inherit scale_adapter return
#'
#' @export
#'
#' @examples
#' proposal <- barker_proposal()
#' adapter <- robust_shape_adapter(initial_scale = 1., target_accept_prob = 0.4)
#' adapter$initialize(proposal, chain_state(c(0, 0)))
robust_shape_adapter <- function(
  initial_scale = NULL, target_accept_prob = NULL, kappa = 0.6
) {
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

is_adapter <- function(object) {
  is.list(object) &&
    all(c("initialize", "update", "finalize", "state") %in% names(object))
}

#' Create a progressive adaptation schedule for use with [sample_chain()].
#'
#' Returns a function (a schedule constructor) that, given the total number of
#' warm-up iterations, produces a three-stage adaptation schedule:
#'
#' * **Stage 1** (`n_fixed_shape_iteration` iterations): adapt scale only,
#'   keeping the proposal shape fixed at the identity. This lets the step size
#'   stabilise before any covariance estimation begins.
#' * **Stage 2** (`n_diagonal_shape_iteration` iterations): adapt scale and
#'   learn a diagonal proposal shape (per-dimension variances).
#' * **Stage 3** (remaining iterations): adapt scale and learn a dense
#'   proposal shape (full covariance matrix).
#'
#' If the sum of `n_fixed_shape_iteration` and `n_diagonal_shape_iteration`
#' exceeds the total number of warm-up iterations, the two counts are reduced
#' proportionally so that each gets half of the available iterations and Stage 3
#' is skipped.
#'
#' @param n_fixed_shape_iteration Number of iterations in Stage 1 (scale only).
#'   Default 50.
#' @param n_diagonal_shape_iteration Number of iterations in Stage 2 (scale +
#'   diagonal shape). Default 50.
#' @param scale_adapter_ Adapter object for the scale. Defaults to
#'   [scale_adapter()].
#' @param diagonal_shape_adapter_ Adapter object for the diagonal shape stage.
#'   Defaults to `shape_adapter("variance")`.
#' @param dense_shape_adapter_ Adapter object for the dense shape stage.
#'   Defaults to `shape_adapter("covariance")`.
#'
#' @return A function that accepts `n_warm_up_iteration` and returns a staged
#'   adapter list suitable for the `adapters` argument of [sample_chain()].
#'
#' @export
#'
#' @examples
#' target_distribution <- list(
#'   log_density = function(x) -sum(x^2) / 2,
#'   gradient_log_density = function(x) -x
#' )
#' withr::with_seed(876287L, {
#'   results <- sample_chain(
#'     target_distribution,
#'     initial_state = stats::rnorm(2),
#'     n_warm_up_iteration = 1000,
#'     n_main_iteration = 1000,
#'     adapters = progressive_adaptation_schedule()
#'   )
#' })
progressive_adaptation_schedule <- function(
  n_fixed_shape_iteration = 50L, # L enforces integer rather than floating-point
  n_diagonal_shape_iteration = 50L,
  scale_adapter_ = scale_adapter(),
  diagonal_shape_adapter_ = shape_adapter("variance"),
  dense_shape_adapter_ = shape_adapter("covariance")
) {
  function(n_warm_up_iteration) {
    # If the two early stages together exceed available iterations, share them
    # equally and skip the dense stage entirely.
    if (n_fixed_shape_iteration + n_diagonal_shape_iteration > n_warm_up_iteration) {
      n_fixed_shape_iteration <- n_warm_up_iteration %/% 2L
      n_diagonal_shape_iteration <- n_warm_up_iteration - n_fixed_shape_iteration
    }
    n_dense_shape_iteration <- (
      n_warm_up_iteration - n_fixed_shape_iteration - n_diagonal_shape_iteration
    )
    stages <- list(
      list(scale_adapter_, n_fixed_shape_iteration),
      list(scale_adapter_, diagonal_shape_adapter_, n_diagonal_shape_iteration)
    )
    if (n_dense_shape_iteration > 0L) {
      stages <- c(stages, list(list(scale_adapter_, dense_shape_adapter_)))
    }
    stages
  }
}
