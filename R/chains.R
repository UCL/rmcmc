#' Sample a Markov chain
#'
#' Sample a Markov chain using Metropolis-Hastings kernel with a user-specified
#' target distribution and proposal (defaulting to Barker proposal), optionally
#' adapting proposal parameters in a warm-up stage.
#'
#' @param target_distribution Target stationary distribution for chain. One of:
#'   * A one-sided formula specifying expression for log density of target
#'     distribution which will be passed to
#'     [target_distribution_from_log_density_formula()] to construct functions
#'     to evaluate log density and its gradient using [deriv()].
#'   * A `bridgestan::StanModel` instance (requires `bridgestan` to be
#'     installed) specifying target model and data. Will be passed to
#'     [target_distribution_from_stan_model()] using default values for optional
#'     arguments - to override call [target_distribution_from_stan_model()]
#'     directly and pass the returned list as the `target_distribution` argument
#'     here.
#'   * A list with named entries `log_density` and `gradient_log_density`
#'     corresponding to respectively functions for evaluating the logarithm of
#'     the (potentially unnormalized) density of the target distribution and its
#'     gradient (only required for gradient-based proposals). As an alternative
#'     to `gradient_log_density` an entry `value_and_gradient_log_density` may
#'     instead be provided which is a function returning both the value and
#'     gradient of the logarithm of the (unnormalized) density of the target
#'     distribution as a list under the names `value` and `gradient`
#'     respectively. The list may also contain a named entry `trace_function`,
#'     correspond to a function which given current chain state outputs a named
#'     vector or list of variables to trace on each main (non-adaptive) chain
#'     iteration. If a `trace_function` entry is not specified, then the default
#'     behaviour is to trace the position component of the chain state along
#'     with the log density of the target distribution.
#' @param initial_state Initial chain state. Either a vector specifying just
#'   the position component of the chain state or a list output by `chain_state`
#'   specifying the full chain state.
#' @param n_warm_up_iteration Number of warm-up (adaptive) chain iterations to
#'   run.
#' @param n_main_iteration Number of main (non-adaptive) chain iterations to
#'   run.
#' @param proposal Proposal distribution object. Defaults to Barker proposal,
#'   that is the output of [barker_proposal()]. Proposal objects are lists which
#'   must minimally define entries `sample`, a function to generate sample from
#'   proposal distribution given current chain state and `log_density_ratio`, a
#'   function to compute log density ratio for proposal for a given pair of
#'   current and proposed chain states. If adapters are being used to adaptively
#'   tune the proposal scale and shape parameters, which is the default
#'   behaviour of `sample_chain`, then additionally the list must also define
#'   entries: `update` a function for updating parameters of proposal,
#'   `parameters` a function for getting current proposal parameter values,
#'   `default_target_accept_prob` a function for getting proposal specific
#'   default target acceptance probability for scale adaptation and
#'   `default_initial_scale` a function for getting proposal and dimension
#'   dependent default initial value for scale parameter.
#' @param adapters List of adapters to tune proposal parameters during warm-up.
#'   Defaults to using list with instances of [scale_adapter()] and
#'   [shape_adapter()], corresponding to respectively, adapting the scale to
#'   coerce the average acceptance rate to a target value using a dual-averaging
#'   algorithm, and adapting the shape to an estimate of the covariance of the
#'   target distribution.
#' @param show_progress_bar Whether to show progress bars during sampling.
#'   Requires `progress` package to be installed to have an effect.
#' @param trace_warm_up Whether to record chain traces and adaptation /
#'   transition statistics during (adaptive) warm-up iterations in addition to
#'   (non-adaptive) main chain iterations.
#'
#' @return A list with entries
#' * `final_state`: the final chain state,
#' * `traces`: a matrix with named columns contained traced variables for each
#'   main chain iteration, with variables along columns and iterations along
#'   rows.
#' * `statistics`: a matrix with named columns containing transition statistics
#'   for each main chain iteration, with statistics along columns and iterations
#'   along rows.
#' * `warm_up_traces`: a matrix with named columns contained traced variables
#'   for each warm-up chain iteration, with variables along columns and
#'   iterations along rows. Only present if `trace_warm_up = TRUE`.
#' * `warm_up_statistics`: a matrix with named columns containing adaptation and
#'   transition statistics for each warm-up chain iteration, with statistics
#'   along columns and iterations along rows. Only present if
#'   `trace_warm_up = TRUE`.
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
#'     n_main_iteration = 1000
#'   )
#' })
sample_chain <- function(
    target_distribution,
    initial_state,
    n_warm_up_iteration,
    n_main_iteration,
    proposal = barker_proposal(),
    adapters = list(scale_adapter(), shape_adapter()),
    show_progress_bar = TRUE,
    trace_warm_up = FALSE) {
  progress_available <- requireNamespace("progress", quietly = TRUE)
  use_progress_bar <- progress_available && show_progress_bar
  initial_state <- check_and_process_initial_state(initial_state)
  target_distribution <- check_and_process_target_distribution(
    target_distribution
  )
  trace_function <- get_trace_function(target_distribution)
  statistic_names <- list("accept_prob")
  warm_up_results <- chain_loop(
    stage_name = "Warm-up",
    n_iteration = n_warm_up_iteration,
    state = initial_state,
    target_distribution = target_distribution,
    proposal = proposal,
    adapters = adapters,
    use_progress_bar = use_progress_bar,
    record_traces_and_statistics = trace_warm_up,
    trace_function = trace_function,
    statistic_names = statistic_names
  )
  main_results <- chain_loop(
    stage_name = "Main",
    n_iteration = n_main_iteration,
    state = warm_up_results$final_state,
    target_distribution = target_distribution,
    proposal = proposal,
    adapters = NULL,
    use_progress_bar = use_progress_bar,
    record_traces_and_statistics = TRUE,
    trace_function = trace_function,
    statistic_names = statistic_names
  )
  if (trace_warm_up) {
    return(combine_stage_results(warm_up_results, main_results))
  } else {
    return(main_results)
  }
}

check_and_process_initial_state <- function(initial_state) {
  if (is.vector(initial_state) && is.atomic(initial_state)) {
    chain_state(initial_state)
  } else if (is.vector(initial_state) && "position" %in% names(initial_state)) {
    initial_state$copy()
  } else {
    stop("initial_state must be a vector or list with an entry named position.")
  }
}

check_and_process_target_distribution <- function(target_distribution) {
  if (inherits(target_distribution, "formula")) {
    target_distribution_from_log_density_formula(target_distribution)
  } else if (inherits(target_distribution, "StanModel")) {
    target_distribution_from_stan_model(target_distribution)
  } else if (
    !is.list(target_distribution) ||
      !("log_density" %in% names(target_distribution))
  ) {
    stop("target_distribution invalid - see documentation for allowable types.")
  } else {
    target_distribution
  }
}

get_trace_function <- function(target_distribution) {
  if (is.null(target_distribution$trace_function)) {
    default_trace_function(target_distribution)
  } else {
    target_distribution$trace_function
  }
}

default_trace_function <- function(target_distribution) {
  function(state) {
    list(
      position = state$position(),
      target_log_density = state$log_density(target_distribution)
    )
  }
}

get_progress_bar <- function(use_progress_bar, n_iteration, label) {
  progress_bar_format <- (
    "%s :percent |:bar| :current/:total [:elapsed<:eta] :tick_rate it/s"
  )
  if (use_progress_bar) {
    progress::progress_bar$new(
      format = sprintf(progress_bar_format, label),
      total = n_iteration,
      clear = FALSE
    )
  } else {
    NULL
  }
}

initialize_traces <- function(trace_function, state, n_iteration) {
  trace_names <- names(unlist(trace_function(state)))
  traces <- matrix(ncol = length(trace_names), nrow = n_iteration)
  colnames(traces) <- trace_names
  traces
}

initialize_statistics <- function(statistic_names, adapters, n_iteration) {
  adapter_statistics <- names(unlist(lapply(adapters, function(a) a$state())))
  statistic_names <- c(statistic_names, adapter_statistics)
  statistics <- matrix(ncol = length(statistic_names), nrow = n_iteration)
  colnames(statistics) <- statistic_names
  statistics
}

initialize_adapters <- function(adapters, proposal, state) {
  for (adapter in adapters) {
    adapter$initialize(proposal, state)
  }
  invisible(adapters)
}

update_adapters <- function(
    adapters, proposal, chain_iteration, state_and_statistics) {
  for (adapter in adapters) {
    adapter$update(proposal, chain_iteration, state_and_statistics)
  }
  invisible(adapters)
}

finalize_adapters <- function(adapters, proposal) {
  for (adapter in adapters) {
    if (!is.null(adapter$finalize)) adapter$finalize(proposal)
  }
  invisible(adapters)
}

chain_loop <- function(
    stage_name,
    n_iteration,
    state,
    target_distribution,
    proposal,
    adapters,
    use_progress_bar,
    record_traces_and_statistics,
    trace_function,
    statistic_names) {
  progress_bar <- get_progress_bar(use_progress_bar, n_iteration, stage_name)
  # Only show 10% increments in progress bar to avoid progress bar updates being
  # a bottleneck when chain iteration rate is high
  tick_amount <- max(n_iteration %/% 10, 1)
  initialize_adapters(adapters, proposal, state)
  if (record_traces_and_statistics) {
    traces <- initialize_traces(trace_function, state, n_iteration)
    statistics <- initialize_statistics(statistic_names, adapters, n_iteration)
  } else {
    traces <- NULL
    statistics <- NULL
  }
  for (chain_iteration in seq_len(n_iteration)) {
    state_and_statistics <- sample_metropolis_hastings(
      state, target_distribution, proposal
    )
    update_adapters(adapters, proposal, chain_iteration, state_and_statistics)
    state <- state_and_statistics$state
    if (record_traces_and_statistics) {
      traces[chain_iteration, ] <- unlist(trace_function(state))
      adapter_states <- lapply(adapters, function(a) a$state())
      statistics[chain_iteration, ] <- unlist(
        c(state_and_statistics$statistics, adapter_states)
      )
    }
    if (!is.null(progress_bar) && (chain_iteration %% tick_amount == 0)) {
      progress_bar$tick(tick_amount)
    }
  }
  # Ensure progress bar shows completed in cases tick_amount not a factor of
  # n_iteration
  if (!is.null(progress_bar) && !progress_bar$finished) progress_bar$update(1)
  finalize_adapters(adapters, proposal)
  list(final_state = state, traces = traces, statistics = statistics)
}

combine_stage_results <- function(warm_up_results, main_results) {
  list(
    final_state = main_results$state,
    traces = main_results$traces,
    statistics = main_results$statistics,
    warm_up_traces = warm_up_results$traces,
    warm_up_statistics = warm_up_results$statistics
  )
}
