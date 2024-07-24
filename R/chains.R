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
  state <- initial_state
  if (is.null(trace_function)) {
    trace_function <- function(state) {
      list(
        position = state$position(),
        target_log_density = state$log_density(target_distribution)
      )
    }
  }
  for (adapter in adapters) {
    adapter$initialize(state)
  }
  progress_bar_format <- (
    "%s :percent |:bar| :current/:total [:elapsed<:eta] :tick_rate it/s")
  if (progress_available && show_progress_bar) {
    progress_bar <- progress::progress_bar$new(
      format = sprintf(progress_bar_format, "Warm-up"),
      total = n_warm_up_iteration,
      clear = FALSE
    )
  }
  for (s in 1:n_warm_up_iteration) {
    state_and_statistics <- sample_metropolis_hastings(
      state, target_distribution, proposal
    )
    for (adapter in adapters) {
      adapter$update(s + 1, state_and_statistics)
    }
    state <- state_and_statistics$state
    if (progress_available && show_progress_bar) progress_bar$tick()
  }
  for (adapter in adapters) {
    if (!is.null(adapter$finalize)) adapter$finalize()
  }
  trace_values <- unlist(trace_function(state))
  traces <- data.frame(matrix(ncol = length(trace_values), nrow = n_main_iteration))
  colnames(traces) <- names(trace_values)
  statistic_values <- unlist(state_and_statistics$statistics)
  statistics <- data.frame(matrix(ncol = length(statistic_values), nrow = n_main_iteration))
  colnames(statistics) <- names(statistic_values)
  if (progress_available && show_progress_bar) {
    progress_bar <- progress::progress_bar$new(
      format = sprintf(progress_bar_format, "Main"),
      total = n_main_iteration,
      clear = FALSE
    )
  }
  for (s in 1:n_main_iteration) {
    state_and_statistics <- sample_metropolis_hastings(
      state, target_distribution, proposal
    )
    state <- state_and_statistics$state
    trace_values <- unlist(trace_function(state))
    traces[s, ] <- trace_values
    statistics[s, ] <- unlist(state_and_statistics$statistics)
    if (progress_available && show_progress_bar) progress_bar$tick()
  }
  list(
    final_state = state,
    traces = traces,
    statistics = statistics
  )
}
