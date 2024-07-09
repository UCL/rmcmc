standard_normal_target_distribution <- function() {
  list(
    log_density = function(x) sum(x^2) / 2,
    grad_log_density = function(x) x
  )
}

default_seed <- function() 9821415L

count_calls <- function(f) {
  call_count <- 0
  wrapped_f <- function(...) {
    call_count <<- call_count + 1
    f(...)
  }
  list(wrapped_f=wrapped_f, get_call_count=function() call_count)
}

check_chain_state <- function(state) {
  expect_type(state, "list")
  expected_names = c(
    "position",
    "momentum",
    "log_density",
    "grad_log_density",
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
  act$n_matching_indices = sum(act$val == dif$val)
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
