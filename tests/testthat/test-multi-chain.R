# Tests for multi-chain related features added for issue #18:
#   * initial_state as a generator function (uses target_distribution$dimension)
#   * chain_index argument and its use in progress messages
#   * combine_chain_results() helper
#   * end-to-end sequential future.apply invocation

# ---------------------------------------------------------------------------
# initial_state as a generator function
# ---------------------------------------------------------------------------

test_that("initial_state as a function uses target_distribution$dimension", {
  dimension <- 3
  target_distribution <- standard_normal_target_distribution(
    dimension = dimension
  )
  adapters <- list(scale_adapter("stochastic_approximation", initial_scale = 1.))
  withr::with_seed(default_seed(), {
    results <- sample_chain(
      target_distribution = target_distribution,
      initial_state = stats::rnorm,
      n_warm_up_iteration = 5,
      n_main_iteration = 10,
      adapters = adapters,
      show_progress_bar = FALSE
    )
  })
  expect_named(
    results,
    c("final_state", "traces", "statistics"),
    ignore.order = TRUE
  )
  # trace matrix should have dimension + 1 columns (positions + log_density)
  expect_ncol(results$traces, dimension + 1)
  expect_nrow(results$traces, 10)
})

test_that("initial_state as a function errors when dimension is missing", {
  target_distribution <- standard_normal_target_distribution()  # no dimension
  expect_error(
    sample_chain(
      target_distribution = target_distribution,
      initial_state = stats::rnorm,
      n_warm_up_iteration = 1,
      n_main_iteration = 1,
      show_progress_bar = FALSE
    ),
    "dimension"
  )
})

test_that("initial_state function returning wrong length is rejected", {
  target_distribution <- standard_normal_target_distribution(dimension = 3)
  bad_generator <- function(n) rep(0, n + 1)  # returns wrong length
  expect_error(
    sample_chain(
      target_distribution = target_distribution,
      initial_state = bad_generator,
      n_warm_up_iteration = 1,
      n_main_iteration = 1,
      show_progress_bar = FALSE
    ),
    "atomic vector of length"
  )
})

test_that("initial_state = runif works and reproduces via with_seed", {
  target_distribution <- standard_normal_target_distribution(dimension = 2)
  run_once <- function() {
    withr::with_seed(42L, {
      sample_chain(
        target_distribution = target_distribution,
        initial_state = stats::runif,
        n_warm_up_iteration = 5,
        n_main_iteration = 5,
        show_progress_bar = FALSE
      )
    })
  }
  r1 <- run_once()
  r2 <- run_once()
  expect_equal(r1$traces, r2$traces)
})

# ---------------------------------------------------------------------------
# chain_index argument and progress-message prefixing
# ---------------------------------------------------------------------------

test_that("chain_index is prefixed onto fallback progress messages", {
  target_distribution <- standard_normal_target_distribution(dimension = 2)
  adapters <- list(scale_adapter("stochastic_approximation", initial_scale = 1.))
  with_mocked_bindings(
    is_progress_package_available = function() FALSE,
    is_progressr_package_available = function() FALSE,
    .package = "rmcmc",
    {
      msgs <- capture_messages(
        sample_chain(
          chain_index = 7L,
          target_distribution = target_distribution,
          initial_state = stats::rnorm,
          n_warm_up_iteration = 10,
          n_main_iteration = 10,
          adapters = adapters,
          show_progress_bar = TRUE
        )
      )
    }
  )
  # every progress line should carry the "Chain 7 |" prefix
  progress_msgs <- grep("done \\(", msgs, value = TRUE)
  expect_true(length(progress_msgs) > 0)
  expect_true(all(grepl("^Chain 7 \\| ", progress_msgs)))
})

test_that("chain_index = NULL leaves progress messages unprefixed", {
  target_distribution <- standard_normal_target_distribution(dimension = 2)
  adapters <- list(scale_adapter("stochastic_approximation", initial_scale = 1.))
  with_mocked_bindings(
    is_progress_package_available = function() FALSE,
    is_progressr_package_available = function() FALSE,
    .package = "rmcmc",
    {
      msgs <- capture_messages(
        sample_chain(
          target_distribution = target_distribution,
          initial_state = stats::rnorm,
          n_warm_up_iteration = 10,
          n_main_iteration = 10,
          adapters = adapters,
          show_progress_bar = TRUE
        )
      )
    }
  )
  progress_msgs <- grep("done \\(", msgs, value = TRUE)
  expect_true(length(progress_msgs) > 0)
  # None of the messages should start with "Chain "
  expect_false(any(grepl("^Chain ", progress_msgs)))
})

test_that("chain_index doubles as future_lapply iteration variable", {
  # future.apply passes X[[i]] as the first positional argument to FUN.
  # We simulate that pattern with plain lapply here (no parallel dependency).
  target_distribution <- standard_normal_target_distribution(dimension = 2)
  adapters <- list(scale_adapter("stochastic_approximation", initial_scale = 1.))
  results <- lapply(
    1:3,
    sample_chain,
    target_distribution = target_distribution,
    initial_state = stats::rnorm,
    n_warm_up_iteration = 5,
    n_main_iteration = 5,
    adapters = adapters,
    show_progress_bar = FALSE
  )
  expect_length(results, 3)
  for (r in results) {
    expect_named(
      r, c("final_state", "traces", "statistics"), ignore.order = TRUE
    )
    expect_nrow(r$traces, 5)
  }
})

# ---------------------------------------------------------------------------
# combine_chain_results()
# ---------------------------------------------------------------------------

test_that("combine_chain_results inverts per-chain nesting", {
  target_distribution <- standard_normal_target_distribution(dimension = 2)
  adapters <- list(scale_adapter("stochastic_approximation", initial_scale = 1.))
  n_chain <- 4
  per_chain <- lapply(
    seq_len(n_chain),
    sample_chain,
    target_distribution = target_distribution,
    initial_state = stats::rnorm,
    n_warm_up_iteration = 3,
    n_main_iteration = 7,
    adapters = adapters,
    show_progress_bar = FALSE
  )
  combined <- combine_chain_results(per_chain)
  expect_named(
    combined,
    c("final_state", "traces", "statistics"),
    ignore.order = TRUE
  )
  for (field in c("final_state", "traces", "statistics")) {
    expect_length(combined[[field]], n_chain)
  }
  # Each traces entry should be a 7-row matrix matching the original
  for (i in seq_len(n_chain)) {
    expect_identical(combined$traces[[i]], per_chain[[i]]$traces)
    expect_identical(combined$statistics[[i]], per_chain[[i]]$statistics)
    expect_identical(combined$final_state[[i]], per_chain[[i]]$final_state)
  }
})

test_that("combine_chain_results preserves warm_up_* fields when present", {
  target_distribution <- standard_normal_target_distribution(dimension = 2)
  adapters <- list(scale_adapter("stochastic_approximation", initial_scale = 1.))
  per_chain <- lapply(1:2, function(i) {
    sample_chain(
      chain_index = i,
      target_distribution = target_distribution,
      initial_state = stats::rnorm,
      n_warm_up_iteration = 4,
      n_main_iteration = 6,
      adapters = adapters,
      trace_warm_up = TRUE,
      show_progress_bar = FALSE
    )
  })
  combined <- combine_chain_results(per_chain)
  expect_named(
    combined,
    c(
      "final_state", "traces", "statistics",
      "warm_up_traces", "warm_up_statistics"
    ),
    ignore.order = TRUE
  )
  expect_length(combined$warm_up_traces, 2)
  expect_length(combined$warm_up_statistics, 2)
})

test_that("combine_chain_results rejects empty input", {
  expect_error(combine_chain_results(list()), "non-empty")
})

test_that("combine_chain_results rejects inconsistent per-chain field names", {
  a <- list(final_state = 1, traces = matrix(0, 1, 1))
  b <- list(final_state = 1, statistics = matrix(0, 1, 1))
  expect_error(combine_chain_results(list(a, b)), "same names")
})

test_that("combine_chain_results rejects unnamed inner lists", {
  a <- list(1, 2)  # no names
  expect_error(combine_chain_results(list(a)), "named lists")
})

# ---------------------------------------------------------------------------
# End-to-end: future.apply with sequential plan
# ---------------------------------------------------------------------------

test_that("future.apply::future_lapply with sequential plan produces N chains", {
  skip_if_not_installed("future.apply")
  target_distribution <- standard_normal_target_distribution(dimension = 2)
  adapters <- list(scale_adapter("stochastic_approximation", initial_scale = 1.))
  old_plan <- future::plan(future::sequential)
  on.exit(future::plan(old_plan), add = TRUE)
  results <- future.apply::future_lapply(
    1:3,
    sample_chain,
    target_distribution = target_distribution,
    initial_state = stats::rnorm,
    n_warm_up_iteration = 5,
    n_main_iteration = 5,
    adapters = adapters,
    show_progress_bar = FALSE,
    future.seed = default_seed()
  )
  expect_length(results, 3)
  combined <- combine_chain_results(results)
  expect_length(combined$traces, 3)
  # Chains started from different random positions should not be identical
  expect_false(identical(combined$traces[[1]], combined$traces[[2]]))
})

# ---------------------------------------------------------------------------
# progressr integration
# ---------------------------------------------------------------------------

test_that("progressr progressor is invoked when progressr is available", {
  skip_if_not_installed("progressr")
  target_distribution <- standard_normal_target_distribution(dimension = 2)
  adapters <- list(scale_adapter("stochastic_approximation", initial_scale = 1.))
  call_count <- 0
  # Stub progressr::progressor to count how often the returned callable is
  # invoked from within chain_loop.
  fake_progressor_factory <- function(steps = 1, message = NULL, ...) {
    function(...) call_count <<- call_count + 1
  }
  local({
    # Only override if progressr is actually available (skip_if handles rest)
    with_mocked_bindings(
      is_progressr_package_available = function() TRUE,
      .package = "rmcmc",
      {
        # Also patch progressr::progressor via a local mock so we don't need a
        # with_progress() wrapper for this unit test.
        orig_progressor <- progressr::progressor
        unlockBinding("progressor", asNamespace("progressr"))
        assign("progressor", fake_progressor_factory, envir = asNamespace("progressr"))
        on.exit({
          assign("progressor", orig_progressor, envir = asNamespace("progressr"))
          lockBinding("progressor", asNamespace("progressr"))
        }, add = TRUE)
        sample_chain(
          target_distribution = target_distribution,
          initial_state = stats::rnorm,
          n_warm_up_iteration = 10,
          n_main_iteration = 10,
          adapters = adapters,
          show_progress_bar = TRUE
        )
      }
    )
  })
  # Expect at least one call per stage (warm-up + main), and typically ~10.
  expect_true(call_count >= 2)
})

test_that("show_progress_bar = FALSE suppresses progressr calls", {
  skip_if_not_installed("progressr")
  target_distribution <- standard_normal_target_distribution(dimension = 2)
  adapters <- list(scale_adapter("stochastic_approximation", initial_scale = 1.))
  call_count <- 0
  fake_progressor_factory <- function(steps = 1, message = NULL, ...) {
    function(...) call_count <<- call_count + 1
  }
  local({
    orig_progressor <- progressr::progressor
    unlockBinding("progressor", asNamespace("progressr"))
    assign("progressor", fake_progressor_factory, envir = asNamespace("progressr"))
    on.exit({
      assign("progressor", orig_progressor, envir = asNamespace("progressr"))
      lockBinding("progressor", asNamespace("progressr"))
    }, add = TRUE)
    sample_chain(
      target_distribution = target_distribution,
      initial_state = stats::rnorm,
      n_warm_up_iteration = 10,
      n_main_iteration = 10,
      adapters = adapters,
      show_progress_bar = FALSE
    )
  })
  expect_equal(call_count, 0)
})
