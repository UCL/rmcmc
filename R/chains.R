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
#' @param adapters Specifies adapters to tune proposal parameters during
#'   warm-up. One of:
#'   * A flat list of adapter objects (e.g. `list(scale_adapter(),
#'     shape_adapter())`), in which case all adapters are active for the full
#'     warm-up period.
#'   * A list of stage specifications, where each stage is itself a list
#'     containing adapter objects and optionally a final integer element giving
#'     the number of iterations for that stage. The last stage may omit the
#'     iteration count, in which case it runs for the remaining warm-up
#'     iterations. For example:
#'     ```
#'     list(
#'       list(scale_adapter(), 50),
#'       list(scale_adapter(), shape_adapter("variance"), 50),
#'       list(scale_adapter(), shape_adapter("covariance"))
#'     )
#'     ```
#'   * A function that accepts `n_warm_up_iteration` as its sole argument and
#'     returns a staged list as above. This is useful for convenience schedule
#'     constructors such as [progressive_adaptation_schedule()].
#'
#'   Defaults to `list(scale_adapter(), shape_adapter())`, which adapts both
#'   scale and covariance shape for the full warm-up period.
#' @param show_progress_bar Whether to show progress bars during sampling. If the
#'   `progress` package is installed, displays a progress bar; otherwise prints
#'   periodic progress messages to the console.
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
  trace_warm_up = FALSE
) {
  progress_available <- is_progress_package_available()
  if (show_progress_bar && !progress_available) {
    message(
      "progress package is not installed, so will print progress updates below."
    )
  }
  initial_state <- check_and_process_initial_state(initial_state)
  target_distribution <- check_and_process_target_distribution(
    target_distribution
  )
  stages <- check_and_process_adapters(adapters, n_warm_up_iteration)
  trace_function <- get_trace_function(target_distribution)
  statistic_names <- list("accept_prob")
  warm_up_results <- list(final_state = initial_state)
  for (stage_index in seq_along(stages)) {
    stage <- stages[[stage_index]]
    stage_warm_up_results <- chain_loop(
      stage_name = sprintf("Warm-up (stage %d/%d)", stage_index, length(stages)),
      n_iteration = stage$n_iteration,
      state = warm_up_results$final_state,
      target_distribution = target_distribution,
      proposal = proposal,
      adapters = stage$adapters,
      show_progress_bar = show_progress_bar,
      progress_available = progress_available,
      record_traces_and_statistics = trace_warm_up,
      trace_function = trace_function,
      statistic_names = statistic_names
    )
    warm_up_results <- combine_warm_up_results(warm_up_results, stage_warm_up_results)
  }
  main_results <- chain_loop(
    stage_name = "Main",
    n_iteration = n_main_iteration,
    state = warm_up_results$final_state,
    target_distribution = target_distribution,
    proposal = proposal,
    adapters = NULL,
    show_progress_bar = show_progress_bar,
    progress_available = progress_available,
    record_traces_and_statistics = TRUE,
    trace_function = trace_function,
    statistic_names = statistic_names
  )
  if (trace_warm_up) {
    combine_stage_results(warm_up_results, main_results)
  } else {
    main_results
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

#' Parse and normalise the adapters argument into a canonical staged list.
#'
#' This internal helper translates flexible user inputs into the strict format
#' required by `chain_loop()`. The `adapters` argument accepts three forms:
#'
#' 1. A flat list of adapter objects -> wraps into a single stage using all
#'    n_warm_up_iteration iterations.
#' 2. A staged list where each element is a list of adapters with an optional
#'    trailing integer giving the stage iteration count. The last stage may
#'    omit the count and will receive the remaining warm-up iterations.
#' 3. A function that accepts n_warm_up_iteration and returns form (2).
#'
#' Returns a list of stages, each a named list with entries:
#'   $adapters    - list of adapter objects for that stage
#'   $n_iteration - integer number of iterations for that stage
#'
#' @noRd
check_and_process_adapters <- function(adapters, n_warm_up_iteration) {
  error_message <- paste(
    "adapters invalid - must be a flat list of adapter objects, a list of",
    "stage specifications (each a list of adapters with an optional trailing",
    "integer iteration count), or a function accepting n_warm_up_iteration."
  )
  if (is.function(adapters)) {
    # Form 3: schedule constructor function
    stages <- adapters(n_warm_up_iteration)
    return(check_and_process_adapters(stages, n_warm_up_iteration))
  }
  if (!is.list(adapters)) stop(error_message)
  # Form 1: flat list of adapter objects
  if (all(sapply(adapters, is_adapter))) {
    return(list(list(adapters = adapters, n_iteration = n_warm_up_iteration)))
  }
  # Form 2: list of stage specifications
  if (!all(sapply(adapters, is.list))) stop(error_message)
  stages <- vector("list", length(adapters))
  n_allocated <- 0L
  for (i in seq_along(adapters)) {
    stage_spec <- adapters[[i]]
    # The last element of a stage spec may be an integer giving n_iteration
    last <- stage_spec[[length(stage_spec)]]
    if (is.numeric(last) && length(last) == 1L && !is_adapter(last)) {
      stage_n_iter <- as.integer(last)
      stage_adapters <- stage_spec[-length(stage_spec)]
    } else {
      # No iteration count supplied — only valid for the last stage
      if (i < length(adapters)) {
        stop(paste(
          "Only the last stage may omit an iteration count.",
          sprintf("Stage %d has no iteration count.", i)
        ))
      }
      stage_n_iter <- n_warm_up_iteration - n_allocated
      if (state_n_iter < 0) stop("Per-stage iteration counts exceeds n_warm_up_iteration")
      stage_adapters <- stage_spec
    }
    if (!all(sapply(stage_adapters, is_adapter))) stop(error_message)
    n_allocated <- n_allocated + stage_n_iter
    stages[[i]] <- list(adapters = stage_adapters, n_iteration = stage_n_iter)
  }
  if (n_allocated != n_warm_up_iteration) {
    stop(sprintf(
      paste(
        "Sum of stage iteration counts (%d) does not equal",
        "n_warm_up_iteration (%d)."
      ),
      n_allocated, n_warm_up_iteration
    ))
  }
  stages
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

is_package_available <- function(pkg) requireNamespace(pkg, quietly = TRUE)

is_progress_package_available <- function() is_package_available("progress")

print_fallback_progress <- function(
  stage_name, chain_iteration, n_iteration, start_time
) {
  elapsed <- proc.time()[["elapsed"]] - start_time
  pct <- round(100 * chain_iteration / n_iteration)
  message(sprintf(
    "%s: %d%% done (%d/%d iterations) | elapsed: %.1fs",
    stage_name, pct, chain_iteration, n_iteration, elapsed
  ))
}

get_progress_bar <- function(
  show_progress_bar, progress_available, n_iteration, label
) {
  progress_bar_format <- (
    "%s :percent |:bar| :current/:total [:elapsed<:eta] :tick_rate it/s"
  )
  if (show_progress_bar && progress_available) {
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
  adapters, proposal, chain_iteration, state_and_statistics
) {
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

chain_loop <- function(  # nolint: cyclocomp_linter. styler: off
  stage_name,
  n_iteration,
  state,
  target_distribution,
  proposal,
  adapters,
  show_progress_bar,
  progress_available,
  record_traces_and_statistics,
  trace_function,
  statistic_names
) {
  progress_bar <- get_progress_bar(
    show_progress_bar, progress_available, n_iteration, stage_name
  )
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
  start_time <- proc.time()[["elapsed"]]
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
    do_tick <- chain_iteration %% tick_amount == 0
    if (!is.null(progress_bar) && do_tick) {
      progress_bar$tick(tick_amount)
    } else if (show_progress_bar && do_tick) { # fallback progress updates
      print_fallback_progress(
        stage_name, chain_iteration, n_iteration, start_time
      )
    }
  }
  # Ensure progress bar shows completed in cases tick_amount not a factor of
  # n_iteration
  progress_unfinished <- n_iteration > 0 && (n_iteration %% tick_amount != 0)
  if (!is.null(progress_bar) && progress_unfinished) {
    progress_bar$update(1)
  } else if (show_progress_bar && progress_unfinished) {
    print_fallback_progress(stage_name, n_iteration, n_iteration, start_time)
  }
  finalize_adapters(adapters, proposal)
  list(final_state = state, traces = traces, statistics = statistics)
}

combine_warm_up_results <- function(warm_up_results_1, warm_up_results_2) {
  # warm_up_results_1$traces and $statistics may be NULL either because this is
  # the first call (initial list(final_state = initial_state)) or because
  # trace_warm_up = FALSE. rbind handles both by treating NULL as empty.
  # 'row bind' stacks matrices on top of each other.
  # rbind(NULL, matrix) just returns the matrix unchanged.
  list(
    final_state = warm_up_results_2$final_state,
    traces = rbind(warm_up_results_1$traces, warm_up_results_2$traces),
    statistics = rbind(warm_up_results_1$statistics, warm_up_results_2$statistics)
  )
}

combine_stage_results <- function(warm_up_results, main_results) {
  list(
    final_state = main_results$final_state,
    traces = main_results$traces,
    statistics = main_results$statistics,
    warm_up_traces = warm_up_results$traces,
    warm_up_statistics = warm_up_results$statistics
  )
}
