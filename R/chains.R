#' Sample a Markov chain
#'
#' Sample a Markov chain using Metropolis-Hastings kernel with given proposal
#' and target distributions, optionally adapting proposal parameters in warm-up
#' stage.
#'
#' @inheritParams sample_metropolis_hastings
#' @param initial_state Initial chain state. Either a vector specifying just
#'   the position component of the chain state or a list output by `chain_state`
#'   specifying the full chain state.
#' @param n_warm_up_iteration Number of warm-up (adaptive) chain iterations to
#'   run.
#' @param n_main_iteration Number of main (non-adaptive) chain iterations to
#'   run.
#' @param adapters List of adapters to tune proposal parameters during warm-up.
#' @param trace_function Function which given current chain state outputs list
#'   of variables to trace on each main (non-adaptive) chain iteration.
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
#' proposal <- barker_proposal(target_distribution, scale = 1.)
#' n_warm_up_iteration <- 1000
#' n_main_iteration <- 1000
#' withr::with_seed(876287L, {
#'   initial_state <- chain_state(stats::rnorm(2))
#'   results <- sample_chain(
#'     target_distribution,
#'     proposal,
#'     initial_state,
#'     n_warm_up_iteration,
#'     n_main_iteration
#'   )
#' })
sample_chain <- function(
    target_distribution,
    proposal,
    initial_state,
    n_warm_up_iteration,
    n_main_iteration,
    adapters = NULL,
    trace_function = NULL,
    show_progress_bar = TRUE,
    trace_warm_up = FALSE) {
  progress_available <- requireNamespace("progress", quietly = TRUE)
  use_progress_bar <- progress_available && show_progress_bar
  if (is.vector(initial_state) && is.atomic(initial_state)) {
    state <- chain_state(initial_state)
  } else if (is.vector(initial_state) && "position" %in% names(initial_state)) {
    state <- initial_state
  } else {
    stop("initial_state must be a vector or list with an entry named position.")
  }
  if (is.null(trace_function)) {
    trace_function <- default_trace_function(target_distribution)
  }
  statistic_names <- list("accept_prob")
  warm_up_results <- chain_loop(
    stage_name = "Warm-up",
    n_iteration = n_warm_up_iteration,
    state = state,
    target_distribution = target_distribution,
    proposal = proposal,
    adapters = adapters,
    use_progress_bar = use_progress_bar,
    record_traces_and_statistics = trace_warm_up,
    trace_function = trace_function,
    statistic_names = statistic_names
  )
  state <- warm_up_results$final_state
  main_results <- chain_loop(
    stage_name = "Main",
    n_iteration = n_main_iteration,
    state = state,
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
    return(
      progress::progress_bar$new(
        format = sprintf(progress_bar_format, label),
        total = n_iteration,
        clear = FALSE
      )
    )
  } else {
    return(NULL)
  }
}

initialize_traces <- function(trace_names, n_iteration) {
  traces <- matrix(ncol = length(trace_names), nrow = n_iteration)
  colnames(traces) <- trace_names
  traces
}

initialize_statistics <- function(statistic_names, n_iteration) {
  statistics <- matrix(ncol = length(statistic_names), nrow = n_iteration)
  colnames(statistics) <- statistic_names
  statistics
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
  for (adapter in adapters) {
    adapter$initialize(state)
  }
  if (record_traces_and_statistics) {
    trace_names <- names(unlist(trace_function(state)))
    traces <- initialize_traces(trace_names, n_iteration)
    adapter_statistics <- names(unlist(lapply(adapters, function(a) a$state())))
    statistics <- initialize_statistics(
      c(statistic_names, adapter_statistics), n_iteration
    )
  } else {
    traces <- NULL
    statistics <- NULL
  }
  for (s in 1:n_iteration) {
    state_and_statistics <- sample_metropolis_hastings(
      state, target_distribution, proposal
    )
    for (adapter in adapters) {
      adapter$update(s + 1, state_and_statistics)
    }
    state <- state_and_statistics$state
    if (record_traces_and_statistics) {
      traces[s, ] <- unlist(trace_function(state))
      adapter_states <- lapply(adapters, function(a) a$state())
      statistics[s, ] <- unlist(
        c(state_and_statistics$statistics, adapter_states)
      )
    }
    if (!is.null(progress_bar)) progress_bar$tick()
  }
  for (adapter in adapters) {
    if (!is.null(adapter$finalize)) adapter$finalize()
  }
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
