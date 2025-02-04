standard_normal_target_distribution <- function(value_and_gradient = FALSE) {
  target_distribution <- list(
    log_density = function(x) -sum(x^2) / 2,
    gradient_log_density = function(x) -x
  )
  if (value_and_gradient) {
    target_distribution[["value_and_gradient_log_density"]] <- function(x) {
      list(
        value = -sum(x^2) / 2, gradient = -x
      )
    }
  }
  target_distribution
}

default_seed <- function() 9821415L

count_calls <- function(f) {
  call_count <- 0
  wrapped_f <- function(...) {
    call_count <<- call_count + 1
    f(...)
  }
  list(wrapped_f = wrapped_f, get_call_count = function() call_count)
}

numerical_gradient <- function(f, h = 1e-8) {
  function(x) {
    apply(
      diag(length(x)), 1, function(e) (f(x + h * e) - f(x - h * e)) / (2 * h)
    )
  }
}

check_chain_state <- function(state) {
  expect_type(state, "list")
  expected_names <- c(
    "position",
    "momentum",
    "log_density",
    "gradient_log_density",
    "dimension",
    "update",
    "copy"
  )
  expect_named(
    state,
    expected_names,
    ignore.order = TRUE,
  )
  for (name in expected_names) {
    expect_type(state[[name]], "closure")
  }
  expect_equal(state$dimension(), length(state$position()))
}

expect_all_different <- function(object, different) {
  act <- quasi_label(rlang::enquo(object), arg = "object")
  dif <- quasi_label(rlang::enquo(different), arg = "different")
  act$n_matching_indices <- sum(act$val == dif$val)
  expect(
    act$n_matching_indices == 0,
    sprintf(
      "%s and %s do not differ in %i indices.",
      act$lab,
      dif$lab,
      act$n_matching_indices
    )
  )
  invisible(act$val)
}

expect_nrow <- function(object, n) {
  act <- quasi_label(rlang::enquo(object), arg = "object")
  act$nrow <- nrow(object)
  expect(
    act$nrow == n,
    sprintf("%s has %i rows not %i.", act$lab, act$nrow, n)
  )
  invisible(act$val)
}

expect_ncol <- function(object, n) {
  act <- quasi_label(rlang::enquo(object), arg = "object")
  act$ncol <- ncol(object)
  expect(
    act$ncol == n,
    sprintf("%s has %i columns not %i.", act$lab, act$ncol, n)
  )
  invisible(act$val)
}

random_covariance_matrix <- function(dimension) {
  temp <- matrix(rnorm(dimension^2), dimension, dimension)
  temp %*% t(temp)
}

multivariate_normal_target_distribution <- function(mean, covariance) {
  chol_covariance <- t(chol(covariance))
  list(
    log_density = function(x) {
      -sum(
        forwardsolve(chol_covariance, x - mean)^2
      ) / 2
    },
    gradient_log_density = function(x) {
      -backsolve(
        t(chol_covariance), forwardsolve(chol_covariance, x - mean)
      )
    }
  )
}

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
