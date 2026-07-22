for (n_warm_up_iteration in c(0, 1, 10)) {
  for (n_main_iteration in c(0, 1, 10, 21)) {
    for (dimension in c(1, 2)) {
      for (trace_warm_up in c(TRUE, FALSE)) {
        for (show_progress_bar in c(TRUE, FALSE)) {
          for (wrapped_initial_state in c(TRUE, FALSE)) {
            for (explicit_trace_function in c(TRUE, FALSE)) {
              test_that(
                sprintf(
                  paste0(
                    "Sampling chain with %i warm-up iterations, ",
                    "%i main iterations, dimension %i, ",
                    "trace_warm_up = %i, show_progress_bar = %i ",
                    "wrapped_initial_state = %i and ",
                    "explicit_trace_function = %i works"
                  ),
                  n_warm_up_iteration,
                  n_main_iteration,
                  dimension,
                  trace_warm_up,
                  show_progress_bar,
                  wrapped_initial_state,
                  explicit_trace_function
                ),
                {
                  target_distribution <- standard_normal_target_distribution()
                  if (explicit_trace_function) {
                    target_distribution[["trace_function"]] <- (
                      default_trace_function(target_distribution)
                    )
                  }
                  adapters <- list(
                    scale_adapter(
                      "stochastic_approximation",
                      initial_scale = 1.
                    )
                  )
                  withr::with_seed(default_seed(), {
                    position <- rnorm(dimension)
                  })
                  if (wrapped_initial_state) {
                    initial_state <- chain_state(position)
                  } else {
                    initial_state <- position
                  }
                  results <- sample_chain(
                    target_distribution = target_distribution,
                    initial_state = initial_state,
                    n_warm_up_iteration = n_warm_up_iteration,
                    n_main_iteration = n_main_iteration,
                    adapters = adapters,
                    trace_warm_up = trace_warm_up,
                    show_progress_bar = show_progress_bar
                  )
                  expected_results_names <- c(
                    "final_state", "traces", "statistics"
                  )
                  if (trace_warm_up) {
                    expected_results_names <- c(
                      expected_results_names,
                      "warm_up_traces",
                      "warm_up_statistics"
                    )
                  }
                  expect_named(
                    results,
                    expected_results_names,
                    ignore.order = TRUE,
                  )
                  expect_nrow(results$traces, n_main_iteration)
                  expect_ncol(results$traces, dimension + 1)
                  expect_nrow(results$statistics, n_main_iteration)
                  expect_ncol(results$statistics, 1)
                  if (trace_warm_up) {
                    expect_nrow(results$warm_up_traces, n_warm_up_iteration)
                    expect_ncol(results$warm_up_traces, dimension + 1)
                    expect_nrow(results$warm_up_statistics, n_warm_up_iteration)
                    expect_ncol(results$warm_up_statistics, 2)
                  }
                }
              )
            }
          }
        }
      }
    }
  }
}

test_that("Sample chains with invalid initial_state raises error", {
  target_distribution <- standard_normal_target_distribution()
  proposal <- barker_proposal()
  expect_error(
    sample_chain(
      target_distribution = target_distribution,
      proposal = proposal,
      initial_state = list(),
      n_warm_up_iteration = 1,
      n_main_iteration = 1,
    ),
    "initial_state"
  )
})

test_that("Sample chains with invalid target_distribution raises error", {
  expect_error(
    sample_chain(
      target_distribution = list(),
      initial_state = c(0., 0.),
      n_warm_up_iteration = 1,
      n_main_iteration = 1,
    ),
    "target_distribution"
  )
})

# ── check_and_process_adapters unit tests ─────────────────────────────────────
#
# check_and_process_adapters is an internal (non-exported) function, so we test
# it indirectly through sample_chain except where we need fine-grained control
# over the returned structure, in which case we call it directly using :::.

test_that(
  "check_and_process_adapters: flat list of adapters -> single stage covering all iterations",
  {
    adapter <- scale_adapter("stochastic_approximation", initial_scale = 1.)
    n_warm_up <- 42L
    stages <- rmcmc:::check_and_process_adapters(list(adapter), n_warm_up)
    expect_length(stages, 1L)
    expect_equal(stages[[1]]$n_iteration, n_warm_up)
    expect_length(stages[[1]]$adapters, 1L)
  }
)

test_that(
  "check_and_process_adapters: staged list with explicit counts on every stage is valid",
  {
    s1 <- scale_adapter("stochastic_approximation", initial_scale = 1.)
    s2 <- scale_adapter("stochastic_approximation", initial_scale = 1.)
    stages <- rmcmc:::check_and_process_adapters(
      list(list(s1, 30L), list(s2, 70L)),
      100L
    )
    expect_length(stages, 2L)
    expect_equal(stages[[1]]$n_iteration, 30L)
    expect_equal(stages[[2]]$n_iteration, 70L)
  }
)

test_that(
  "check_and_process_adapters: last stage may omit iteration count and receives remainder",
  {
    s1 <- scale_adapter("stochastic_approximation", initial_scale = 1.)
    s2 <- scale_adapter("stochastic_approximation", initial_scale = 1.)
    stages <- rmcmc:::check_and_process_adapters(
      list(list(s1, 30L), list(s2)),
      100L
    )
    expect_length(stages, 2L)
    expect_equal(stages[[1]]$n_iteration, 30L)
    expect_equal(stages[[2]]$n_iteration, 70L)
  }
)

test_that(
  "check_and_process_adapters: non-last stage omitting iteration count raises error",
  {
    s1 <- scale_adapter("stochastic_approximation", initial_scale = 1.)
    s2 <- scale_adapter("stochastic_approximation", initial_scale = 1.)
    expect_error(
      rmcmc:::check_and_process_adapters(
        # Stage 1 has no count but is not the last stage
        list(list(s1), list(s2, 50L)),
        100L
      ),
      "Only the last stage"
    )
  }
)

test_that(
  "check_and_process_adapters: stage counts summing to less than n_warm_up_iteration raises error",
  {
    s1 <- scale_adapter("stochastic_approximation", initial_scale = 1.)
    s2 <- scale_adapter("stochastic_approximation", initial_scale = 1.)
    expect_error(
      rmcmc:::check_and_process_adapters(
        list(list(s1, 30L), list(s2, 40L)),
        100L
      ),
      "does not equal"
    )
  }
)

test_that(
  "check_and_process_adapters: stage counts summing to more than n_warm_up_iteration raises error",
  {
    s1 <- scale_adapter("stochastic_approximation", initial_scale = 1.)
    s2 <- scale_adapter("stochastic_approximation", initial_scale = 1.)
    expect_error(
      rmcmc:::check_and_process_adapters(
        list(list(s1, 60L), list(s2, 60L)),
        100L
      ),
      "does not equal"
    )
  }
)

test_that(
  "check_and_process_adapters: function form is called with n_warm_up_iteration and result is validated",
  {
    s1 <- scale_adapter("stochastic_approximation", initial_scale = 1.)
    schedule_fn <- function(n) list(list(s1, n))
    stages <- rmcmc:::check_and_process_adapters(schedule_fn, 77L)
    expect_length(stages, 1L)
    expect_equal(stages[[1]]$n_iteration, 77L)
  }
)

test_that(
  "check_and_process_adapters: non-list non-function input raises error",
  {
    expect_error(
      rmcmc:::check_and_process_adapters("not_valid", 10L),
      "adapters invalid"
    )
  }
)

test_that(
  "check_and_process_adapters: numeric scalar input raises error",
  {
    # A numeric value is neither a list, a function, nor a flat list of adapters
    # nor a list of stage specifications, so should be rejected with the
    # informative "adapters invalid" error.
    expect_error(
      rmcmc:::check_and_process_adapters(42, 10L),
      "adapters invalid"
    )
  }
)

test_that(
  "check_and_process_adapters: staged list with per-stage counts exceeding n_warm_up_iteration and no trailing count raises informative error",
  {
    s1 <- scale_adapter("stochastic_approximation", initial_scale = 1.)
    s2 <- scale_adapter("stochastic_approximation", initial_scale = 1.)
    # Stages 1 and 2 together already consume more than n_warm_up_iteration, and
    # stage 3 has no explicit count so it tries to compute the remainder, which
    # would be negative. We should get a targeted error message *at parse time*
    # rather than a cryptic failure later when attempting to sample a stage
    # with a negative iteration count.
    expect_error(
      rmcmc:::check_and_process_adapters(
        list(list(s1, 60L), list(s2, 60L), list(s1)),
        100L
      ),
      "Per-stage iteration counts exceeds n_warm_up_iteration"
    )
  }
)

test_that(
  "check_and_process_adapters: staged list with an empty middle stage (no adapters, just an iteration count) is valid",
  {
    # A stage with zero adapters but an explicit iteration count represents a
    # "pause" during warm-up where the chain still moves but no proposal
    # parameters are updated. The parser should accept it and allocate the
    # requested number of iterations to that stage.
    s <- scale_adapter("stochastic_approximation", initial_scale = 1.)
    stages <- rmcmc:::check_and_process_adapters(
      list(list(s, 30L), list(50L), list(s)),
      100L
    )
    expect_length(stages, 3L)
    expect_equal(stages[[1]]$n_iteration, 30L)
    expect_length(stages[[1]]$adapters, 1L)
    expect_equal(stages[[2]]$n_iteration, 50L)
    expect_length(stages[[2]]$adapters, 0L)
    expect_equal(stages[[3]]$n_iteration, 20L)
    expect_length(stages[[3]]$adapters, 1L)
  }
)

test_that(
  "check_and_process_adapters: empty top-level list is treated as a single no-adapter stage covering all warm-up iterations",
  {
    # An empty adapters argument means "run warm-up with no adaptation at all".
    # The parser should produce a single stage with zero adapters and
    # n_iteration equal to n_warm_up_iteration.
    stages <- rmcmc:::check_and_process_adapters(list(), 25L)
    expect_length(stages, 1L)
    expect_equal(stages[[1]]$n_iteration, 25L)
    expect_length(stages[[1]]$adapters, 0L)
  }
)

test_that(
  "check_and_process_adapters: staged list containing invalid element raises error",
  {
    s1 <- scale_adapter("stochastic_approximation", initial_scale = 1.)
    expect_error(
      rmcmc:::check_and_process_adapters(
        # Second stage spec has a string instead of an adapter
        list(list(s1, 50L), list("not_an_adapter", 50L)),
        100L
      ),
      "adapters invalid"
    )
  }
)

# ── sample_chain with staged adapters integration tests ───────────────────────

test_that(
  "sample_chain with an empty middle stage (no adapters, just an iteration count) runs without error",
  {
    # Integration-level check that the whole pipeline (parse -> sample -> combine)
    # handles a stage with zero adapters. The empty stage exercises the
    # `for (adapter in adapters)` no-op path in initialize/update/finalize_adapters
    # and the empty-adapter-set path in initialize_statistics.
    target_distribution <- standard_normal_target_distribution()
    s <- scale_adapter("stochastic_approximation", initial_scale = 1.)
    withr::with_seed(default_seed(), {
      expect_no_error(
        results <- sample_chain(
          target_distribution = target_distribution,
          initial_state = rnorm(2),
          n_warm_up_iteration = 20L,
          n_main_iteration = 5L,
          adapters = list(list(s, 5L), list(10L), list(s)),
          trace_warm_up = FALSE,
          show_progress_bar = FALSE
        )
      )
    })
    expect_nrow(results$traces, 5L)
  }
)

test_that(
  "sample_chain with two-stage adapters (same adapter type) produces correct warm-up row counts",
  {
    # Both stages use only a scale adapter so the statistics columns match
    # across stages, making rbind in combine_warm_up_results safe.
    target_distribution <- standard_normal_target_distribution()
    n_warm_up <- 10L
    n_stage_1 <- 4L
    n_stage_2 <- n_warm_up - n_stage_1
    s1 <- scale_adapter("stochastic_approximation", initial_scale = 1.)
    s2 <- scale_adapter("stochastic_approximation", initial_scale = 1.)
    withr::with_seed(default_seed(), {
      results <- sample_chain(
        target_distribution = target_distribution,
        initial_state = rnorm(2),
        n_warm_up_iteration = n_warm_up,
        n_main_iteration = 5L,
        adapters = list(list(s1, n_stage_1), list(s2)),
        trace_warm_up = TRUE,
        show_progress_bar = FALSE
      )
    })
    expect_nrow(results$warm_up_traces, n_warm_up)
    expect_nrow(results$warm_up_statistics, n_warm_up)
  }
)

test_that(
  "sample_chain with three-stage adapters (same adapter type) produces correct warm-up row counts",
  {
    target_distribution <- standard_normal_target_distribution()
    n_warm_up <- 12L
    s1 <- scale_adapter("stochastic_approximation", initial_scale = 1.)
    s2 <- scale_adapter("stochastic_approximation", initial_scale = 1.)
    s3 <- scale_adapter("stochastic_approximation", initial_scale = 1.)
    withr::with_seed(default_seed(), {
      results <- sample_chain(
        target_distribution = target_distribution,
        initial_state = rnorm(2),
        n_warm_up_iteration = n_warm_up,
        n_main_iteration = 5L,
        adapters = list(list(s1, 3L), list(s2, 4L), list(s3)),
        trace_warm_up = TRUE,
        show_progress_bar = FALSE
      )
    })
    expect_nrow(results$warm_up_traces, n_warm_up)
    expect_nrow(results$warm_up_statistics, n_warm_up)
  }
)

test_that(
  "sample_chain with staged adapters (different adapter types) runs without error when trace_warm_up = FALSE",
  {
    # When stages have different adapter sets the statistics matrices have
    # different column counts. With trace_warm_up = FALSE no rbind of statistics
    # occurs, so this should run cleanly.
    target_distribution <- standard_normal_target_distribution()
    s_scale <- scale_adapter("stochastic_approximation", initial_scale = 1.)
    s_shape <- shape_adapter("variance")
    withr::with_seed(default_seed(), {
      expect_no_error(
        sample_chain(
          target_distribution = target_distribution,
          initial_state = rnorm(2),
          n_warm_up_iteration = 10L,
          n_main_iteration = 5L,
          adapters = list(list(s_scale, 4L), list(s_scale, s_shape)),
          trace_warm_up = FALSE,
          show_progress_bar = FALSE
        )
      )
    })
  }
)

test_that(
  "sample_chain with a function as adapters argument runs without error",
  {
    target_distribution <- standard_normal_target_distribution()
    s <- scale_adapter("stochastic_approximation", initial_scale = 1.)
    schedule_fn <- function(n) list(list(s, n %/% 2L), list(s))
    withr::with_seed(default_seed(), {
      expect_no_error(
        sample_chain(
          target_distribution = target_distribution,
          initial_state = rnorm(2),
          n_warm_up_iteration = 10L,
          n_main_iteration = 5L,
          adapters = schedule_fn,
          trace_warm_up = FALSE,
          show_progress_bar = FALSE
        )
      )
    })
  }
)

test_that(
  "sample_chain with staged adapters returns correct final_state type",
  {
    target_distribution <- standard_normal_target_distribution()
    s1 <- scale_adapter("stochastic_approximation", initial_scale = 1.)
    s2 <- scale_adapter("stochastic_approximation", initial_scale = 1.)
    withr::with_seed(default_seed(), {
      results <- sample_chain(
        target_distribution = target_distribution,
        initial_state = rnorm(2),
        n_warm_up_iteration = 10L,
        n_main_iteration = 5L,
        adapters = list(list(s1, 4L), list(s2)),
        trace_warm_up = FALSE,
        show_progress_bar = FALSE
      )
    })
    # final_state should be a chain state object (list with a $position entry)
    expect_true("position" %in% names(results$final_state))
  }
)

test_that(
  "sample_chain with invalid adapters argument raises error",
  {
    target_distribution <- standard_normal_target_distribution()
    expect_error(
      sample_chain(
        target_distribution = target_distribution,
        initial_state = c(0., 0.),
        n_warm_up_iteration = 10L,
        n_main_iteration = 5L,
        adapters = "not_valid"
      ),
      "adapters invalid"
    )
  }
)
make_fallback_test_inputs <- function() {
  target_distribution <- standard_normal_target_distribution()
  adapters <- list(scale_adapter("stochastic_approximation", initial_scale = 1.))
  withr::with_seed(default_seed(), {
    position <- rnorm(2)
  })
  list(
    target_distribution = target_distribution,
    adapters = adapters,
    position = position
  )
}

test_that("Manual progress fallback prints messages when progress unavailable", {
  inputs <- make_fallback_test_inputs()
  n_warm_up_iteration <- 10
  # use non-multiple of 10 to test finalisation of progress updates
  n_main_iteration <- 21
  # Simulate progress package being unavailable by mocking
  with_mocked_bindings(
    is_progress_package_available = function() FALSE,
    .package = "rmcmc",
    {
      msgs <- capture_messages(
        sample_chain(
          target_distribution = inputs$target_distribution,
          initial_state = inputs$position,
          n_warm_up_iteration = n_warm_up_iteration,
          n_main_iteration = n_main_iteration,
          adapters = inputs$adapters,
          show_progress_bar = TRUE
        )
      )
    }
  )
  expected_n_progress_messages <- function(n_iteration) {
    tick_amount <- max(n_iteration %/% 10, 1)
    n_iteration %/% tick_amount + (n_iteration %% tick_amount != 0)
  }
  # 1 upfront warning + n_iteration dependent number of messages (warm-up + main)
  expect_length(
    msgs,
    1 + expected_n_progress_messages(n_warm_up_iteration)
      + expected_n_progress_messages(n_main_iteration)
  )
  expect_true(any(grepl("progress package is not installed", msgs)))
  expect_true(any(grepl("10%", msgs)))
  expect_true(any(grepl("100%", msgs)))
})

test_that("No manual progress output when show_progress_bar is FALSE", {
  inputs <- make_fallback_test_inputs()
  with_mocked_bindings(
    is_progress_package_available = function() FALSE,
    .package = "rmcmc",
    {
      expect_no_message(
        sample_chain(
          target_distribution = inputs$target_distribution,
          initial_state = inputs$position,
          n_warm_up_iteration = 10,
          n_main_iteration = 10,
          adapters = inputs$adapters,
          show_progress_bar = FALSE
        )
      )
    }
  )
})
