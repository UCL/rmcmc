#' Construct a new chain state.
#'
#' The chain state object provides cached access to target distribution
#' log density and its gradient.
#'
#' @param position Position component of chain state.
#' @param momentum Momentum component of chain state. Optional.
#'
#' @return New chain state object. A list with entries
#' * `position`: A zero-argument function to evaluate position vector.
#' * `momentum`: A zero-argument function to evaluate momentum vector.
#' * `dimension`: A zero-argument function evaluate dimension of position and
#'   momentum vectors.
#' * `update`: A function accepting arguments `position` and `momentum` for
#'   updating the value of one or both of these state components.
#' * `copy`: A function for creating a copy of the state object including any
#'   cached values.
#' * `log_density`: A function accepting argument `target_distribution` for
#'   evaluating the log density of the target distribution at the current
#'   state, with caching of the value to avoid recomputation on subsequent
#'   calls.
#' * `grad_log_density`: A function accepting argument `target_distribution` for
#'   evaluating the gradient of the log density of the target distribution at
#'   the current state, with caching of the value to avoid recomputation on
#'   subsequent calls.
#' @export
#'
#' @examples
#' state <- chain_state(c(0.1, -0.5))
#' target_distribution <- list(
#'   log_density = function(x) sum(x^2) / 2,
#'   grad_log_density = function(x) x
#' )
#' state$grad_log_density(target_distribution)
chain_state <- function(position, momentum = NULL) {
  new_chain_state(position, momentum)
}

new_chain_state <- function(
    position,
    momentum = NULL,
    cached_log_density = NULL,
    cached_grad_log_density = NULL) {
  log_density <- function(target_distribution) {
    if (is.null(cached_log_density)) {
      cached_log_density <<- target_distribution$log_density(position)
    }
    cached_log_density
  }
  grad_log_density <- function(target_distribution) {
    if (is.null(cached_grad_log_density)) {
      if ("value_and_grad_log_density" %in% names(target_distribution)) {
        value_and_grad <- target_distribution$value_and_grad_log_density(position)
        cached_log_density <<- value_and_grad$value
        cached_grad_log_density <<- value_and_grad$grad
      } else {
        cached_grad_log_density <<- target_distribution$grad_log_density(position)
      }
    }
    cached_grad_log_density
  }
  update <- function(position = NULL, momentum = NULL) {
    if (!is.null(position)) {
      position <<- position
      cached_log_density <<- NULL
      cached_grad_log_density <<- NULL
    }
    if (!is.null(momentum)) {
      momentum <<- momentum
    }
  }
  list(
    position = function() position,
    momentum = function() momentum,
    dimension = function() length(position),
    update = update,
    copy = function() {
      new_chain_state(
        position, momentum, cached_log_density, cached_grad_log_density
      )
    },
    log_density = log_density,
    grad_log_density = grad_log_density
  )
}
