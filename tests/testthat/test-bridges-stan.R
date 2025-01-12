skip_on_cran()
skip_on_os("windows")

cached_example_gaussian_stan_model <- (
  function() {
    model <- NULL
    function() {
      if (is.null(model)) {
        model <<- example_gaussian_stan_model()
      }
      model
    }
  })()

test_that("Creating example Gaussian Stan model works", {
  model <- cached_example_gaussian_stan_model()
  expect_s3_class(model, c("StanModel", "R6"))
  expect_identical(model$param_names(), c("mu", "sigma"))
  expect_identical(model$param_unc_num(), 2L)
})

test_that("Creating target distribution from Stan model works", {
  model <- cached_example_gaussian_stan_model()
  target_distribution <- target_distribution_from_stan_model(model)
  check_target_distribution(target_distribution)
  position <- rep(0, model$param_unc_num())
  check_log_density_and_gradient(
    position, target_distribution, model$log_density
  )
})

for (include_log_density in c(TRUE, FALSE)) {
  test_that(
    sprintf(
      "Creating trace function from Stan model works (include_log_density = %i)",
      include_log_density
    ),
    {
      model <- cached_example_gaussian_stan_model()
      trace_function <- target_distribution_from_stan_model(
        model,
        include_log_density = include_log_density
      )$trace_function
      expect_type(trace_function, "closure")
      position <- rep(0, model$param_unc_num())
      state <- chain_state(position)
      trace_values <- trace_function(state)
      if (include_log_density) {
        expected_length <- model$param_num() + 1L
      } else {
        expected_length <- model$param_num()
      }
      expect_identical(length(trace_values), expected_length)
      expected_parameter_values <- model$param_constrain(position)
      expected_parameter_names <- model$param_names()
      if (include_log_density) {
        expected_names <- c(expected_parameter_names, "log_density")
      } else {
        expected_names <- expected_parameter_names
      }
      expect_named(trace_values, expected_names)
      for (i in seq_along(expected_parameter_values)) {
        name <- expected_parameter_names[[i]]
        expect_identical(trace_values[[name]], expected_parameter_values[[i]])
      }
      if (include_log_density) {
        expected_log_density <- model$log_density(position)
        expect_identical(trace_values[["log_density"]], expected_log_density)
      }
    }
  )
}

test_that("Using sample_chains with Stan model as target_distribution works", {
  model <- cached_example_gaussian_stan_model()
  dimension <- model$param_unc_num()
  n_warm_up_iteration <- 50
  n_main_iteration <- 50
  results <- sample_chain(
    target_distribution = model,
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
