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

#' Check whether an object is a non-scalar vector
#'
#' @param obj Object to check.
#'
#' @noRd
#'
#' @return `TRUE` if `obj` is a non-scalar vector and `FALSE` otherwise.
is_non_scalar_vector <- function(obj) {
  is.null(dim(obj)) && length(obj) > 1
}

#' Matrix vector multiplication like operator with vectors and scalars
#' considered as diagonal matrices.
#'
#' At least one of arguments must be a vector.
#'
#' @param left Left operand in multiplication. If a scalar considered to be
#'   equivalent to scaled identity matrix. If a vector considered to be
#'   equivalent to a diagonal matrix with vector values along diagonal.
#' @param right Right operand in multiplication. If a scalar considered to be
#'   equivalent to scaled identity matrix. If a vector considered to be
#'   equivalent to a diagonal matrix with vector values along diagonal.
#'
#' @noRd
#'
#' @return Result of matrix vector multiplication of `left` and `right`.
`%@%` <- function(left, right) {
  if (is.null(dim(left)) && is.null(dim(right))) {
    return(left * right)
  } else if (is.matrix(left) && is_non_scalar_vector(right)) {
    return(Matrix::drop(left %*% right))
  } else if (is_non_scalar_vector(left) && is.matrix(right)) {
    return(Matrix::drop(Matrix::t(right) %*% left))
  } else {
    stop("Expected at least one vector argument")
  }
}
