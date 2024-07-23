#' Logistic sigmoid function `1 / (1 + exp(-x))`.
#'
#' @param x Scalar or vector evaluate logistic sigmoid at.
#'
#' @noRd
#'
#' @return Value of logistic sigmoid at `x`.
logistic_sigmoid <- function(x) {
  stats::plogis(x)
}

#' Numerically stable computation of `log(1 + exp(x))`.
#'
#' @param x Scalar or vector to evaluate function at.
#'
#' @noRd
#'
#' @return Value of `log(1 + exp(x))`
log1p_exp <- function(x) {
  pmax(x, 0) + log1p(exp(-abs(x)))
}

#' Numerically stable computation of `min(1, exp(x))`.
#'
#' @param x Scalar or vector to evaluate function at.
#'
#' @noRd
#'
#' @return Value of `min(1, exp(x))`
min_1_exp <- function(x) {
  exp(pmin(0, x))
}
