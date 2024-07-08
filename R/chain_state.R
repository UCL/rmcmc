#' Construct a new chain state.
#'
#' The chain state object provides cached access to target distribution
#' log density and its gradient.
#'
#' @param position Position component of chain state.
#' @param momentum Momentum component of chain state. Optional.
#'
#' @return New chain state object.
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
      if ("value_and_grad_log_density" %in% target_distribution) {
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
