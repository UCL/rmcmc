#' Sample new state from Barker proposal.
#'
#' @param x Current chain state.
#' @param grad_x Gradient of target distribution log density function at current
#'   chain state.
#' @param runif Function which generates a random vector from standard uniform
#'   distribution given an integer size.
#' @param raux Function which generates a random vector from auxiliary variable
#'   distribution.
#' @param ... Any additional arguments to pass to `raux` function.
#'
#' @return Proposed new state.
#' @export
#'
#' @examples
#' x <- c(0., 0.)
#' grad_x <- c(-0.1, 0.5)
#' sample_barker(x, grad_x, sd = 2.)
sample_barker <- function(
    x,
    grad_x,
    raux = stats::rnorm,
    runif = stats::runif,
    ...) {
  n <- length(grad_x)
  z <- raux(n, ...)
  b <- 2 * (runif(n) < logistic_sigmoid(grad_x * z)) - 1
  x + b * z
}

#' Compute logarithm of Barker proposal density ratio.
#'
#' @param x Current chain state.
#' @param y Proposed chain state.
#' @param grad_x Gradient of target distribution log density function at current
#'   chain state.
#' @param grad_y Gradient of target distribution log density function at
#'   proposed chain state.
#'
#' @return Logarithm of proposal density ratio.
#' @export
#'
#' @examples
#' x <- c(0., 0.)
#' grad_x <- c(-0.1, 0.5)
#' y <- c(1., 0.)
#' grad_y <- c(0.7, -0.4)
#' log_density_ratio_barker(x, y, grad_x, grad_y)
log_density_ratio_barker <- function(x, y, grad_x, grad_y) {
  log1p_exp((x - y) * grad_x) - log1p_exp((y - x) * grad_y)
}
