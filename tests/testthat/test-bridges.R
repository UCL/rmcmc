test_that("Constructing target distribution from log density formula works", {
  log_density_formula <- ~ -(x^2 + y^2 + z^2) / 2
  target_distribution <- target_distribution_from_log_density_formula(
    log_density_formula
  )
  check_target_distribution(target_distribution)
  withr::with_seed(seed = default_seed(), position <- rnorm(3))
  check_log_density_and_gradient(
    position, target_distribution, function(x) -sum(x^2) / 2
  )
  trace_values <- target_distribution$trace_function(chain_state(position))
  expect_named(trace_values, c("x", "y", "z", "log_density"))
  expect_equal(
    trace_values, c(position, -sum(position^2) / 2),
    ignore_attr = TRUE
  )
})

test_that("Using sample_chains with formula as target_distribution works", {
  log_density_formula <- ~ -(x^2 + y^2 + z^2) / 2
  dimension <- 3
  n_warm_up_iteration <- 50
  n_main_iteration <- 50
  results <- sample_chain(
    target_distribution = log_density_formula,
    initial_state = rep(0, dimension),
    n_warm_up_iteration = n_warm_up_iteration,
    n_main_iteration = n_main_iteration
  )
  expect_named(
    results,
    c("final_state", "traces", "statistics"),
    ignore.order = TRUE,
  )
  expect_nrow(results$traces, n_main_iteration)
  expect_ncol(results$traces, dimension + 1)
  expect_nrow(results$statistics, n_main_iteration)
  expect_ncol(results$statistics, 1)
})
