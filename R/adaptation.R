step_size_adapter <- function(
    proposal, initial_step_size, target_accept_prob = 0.4, kappa = 0.6) {
  log_step_size <- NULL
  initialize <- function(initial_state) {
    log_step_size <<- log(initial_step_size)
    proposal$update(step_size=initial_step_size)
  }
  update <- function(sample_index, state_and_statistics) {
    gamma <- sample_index^(-kappa)
    accept_prob <- state_and_statistics$statistics$accept_prob
    log_step_size <<- log_step_size + gamma * (accept_prob - target_accept_prob)
    proposal$update(step_size = exp(log_step_size))
  }
  list(initialize = initialize, update = update, finalize = function () {})
}

dense_covariance_adapter <- function(proposal, kappa = 0.6) {
  mean_estimate <- 0
  covariance_estimate <- 0
  outer_difference_from_mean <- function(position) {
    diff <- position - mean_estimate
    outer(diff, diff)
  }
  initialize <- function(initial_state) {
    mean_estimate <<- initial_state$position()
    diff <- initial_state$position() - mean_estimate
    covariance_estimate <<- outer_difference_from_mean(initial_state$position)
  }
  update <- function(sample_index, state_and_statistics) {
    gamma <- sample_index^(-kappa)
    position <- state_and_statistics$state$position()
    mean_estimate <<- mean_estimate + gamma * (position - mean_estimate)
    covariance_estimate <<- covariance_estimate + gamma * (
      outer_difference_from_mean(position) - covariance_estimate)
  }
  finalize <- function() {
    proposal$update(preconditioner=covariance_estimate)
  }
  list(initialize = initialize, update = update, finalize = finalize)
}
