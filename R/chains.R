#' Sample a Markov chain
#'
#' Sample a Markov chain using Metropolis-Hastings kernel with given proposal
#' and target distributions, optionally adapting proposal parameters in warm-up
#' stage.
#'
#' @inheritParams sample_metropolis_hastings
#' @param initial_state Initial chain state.
#' @param n_warm_up_iteration Number of warm-up (adaptive) chain iterations to
#'   run.
#' @param n_main_iteration Number of main (non-adaptive) chain iterations to
#'   run.
#' @param adapters List of adapters to tune proposal parameters during warm-up.
#' @param trace_function Function which given current chain state outputs list
#'   of variables to trace on each main (non-adaptive) chain iteration.
#' @param show_progress_bar Whether to show progress bars during sampling.
#'   Requires `progress` package to be installed to have an effect.
#'
#' @return A list with entries
#' * `final_state`: the final chain state,
#' * `traces`: a dataframe contained traced variables for each main chain
#'   iteration, with variables along columns and iterations along rows.
#' * `statistics`: a dataframe containing transition statistics for each main
#'   chain iteration, with statistics along columns and iterations along rows.
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
    show_progress_bar = TRUE) {
  progress_available <- requireNamespace("progress", quietly = TRUE)
  use_progress_bar <- progress_available && show_progress_bar
  if (is.null(trace_function)) {
    trace_function <- default_trace_function(target_distribution)
  }
  state_and_statistics <- warm_up_chain_loop(
    n_iteration=n_warm_up_iteration,
    state=initial_state,
    target_distribution=target_distribution,
    proposal=proposal,
    adapters=adapters,
    use_progress_bar=use_progress_bar
  )
  main_chain_loop(
    n_iteration=n_main_iteration,
    state_and_statistics=state_and_statistics,
    target_distribution=target_distribution,
    proposal=proposal,
    trace_function=trace_function,
    use_progress_bar=use_progress_bar
  )
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

initialize_traces <- function(trace_values, n_iteration) {
  traces <- data.frame(
    matrix(ncol = length(trace_values), nrow = n_iteration)
  )
  colnames(traces) <- names(trace_values)
  traces
}

initialize_statistics <- function(statistic_values, n_iteration) {
  statistics <- data.frame(
    matrix(ncol = length(statistic_values), nrow = n_iteration)
  )
  colnames(statistics) <- names(statistic_values)
  statistics
}

warm_up_chain_loop <- function(
    n_iteration,
    state,
    target_distribution,
    proposal,
    adapters,
    use_progress_bar) {
  progress_bar <- get_progress_bar(
    use_progress_bar, n_iteration, "Warm-up"
  )
  for (adapter in adapters) {
    adapter$initialize(state)
  }
  for (s in 1:n_iteration) {
    state_and_statistics <- sample_metropolis_hastings(
      state, target_distribution, proposal
    )
    for (adapter in adapters) {
      adapter$update(s + 1, state_and_statistics)
    }
    state <- state_and_statistics$state
    if (!is.null(progress_bar)) progress_bar$tick()
  }
  for (adapter in adapters) {
    if (!is.null(adapter$finalize)) adapter$finalize()
  }
  state_and_statistics
}

main_chain_loop <- function(
    n_iteration,
    state_and_statistics,
    target_distribution,
    proposal,
    trace_function,
    use_progress_bar) {
  state <- state_and_statistics$state
  statistic_values <- unlist(state_and_statistics$statistics)
  trace_values <- unlist(trace_function(state))
  traces <- initialize_traces(trace_values, n_iteration)
  statistics <- initialize_statistics(statistic_values, n_iteration)
  progress_bar <- get_progress_bar(use_progress_bar, n_iteration, "Main")
  for (s in 1:n_iteration) {
    state_and_statistics <- sample_metropolis_hastings(
      state, target_distribution, proposal
    )
    state <- state_and_statistics$state
    trace_values <- unlist(trace_function(state))
    traces[s, ] <- trace_values
    statistics[s, ] <- unlist(state_and_statistics$statistics)
    if (use_progress_bar) progress_bar$tick()
  }
  list(
    final_state = state_and_statistics$state,
    traces = traces,
    statistics = statistics
  )
}
