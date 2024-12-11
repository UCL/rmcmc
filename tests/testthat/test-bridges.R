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

check_target_distribution <- function(target_distribution) {
  expect_type(target_distribution, "list")
  expect_named(
    target_distribution,
    c("log_density", "value_and_gradient_log_density", "trace_function")
  )
  expect_type(target_distribution$log_density, "closure")
  expect_type(target_distribution$value_and_gradient_log_density, "closure")
  expect_type(target_distribution$trace_function, "closure")
}

check_log_density_and_gradient <- function(
    position, target_distribution, true_log_density) {
  log_density <- target_distribution$log_density(position)
  value_and_gradient_log_density <- (
    target_distribution$value_and_gradient_log_density(position)
  )
  expect_type(log_density, "double")
  expect_equal(log_density, true_log_density(position))
  expect_type(value_and_gradient_log_density, "list")
  expect_identical(value_and_gradient_log_density$value, log_density)
  expect_equal(
    value_and_gradient_log_density$gradient,
    numerical_gradient(true_log_density)(position),
    tolerance = 1e-6,
    ignore_attr = TRUE
  )
}

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
