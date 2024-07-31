is_non_scalar_vector <- function(obj) {
  is.null(dim(obj)) && length(obj) > 1
}

get_shape_matrix <- function(scale, shape) {
  shape <- if (is_non_scalar_vector(shape)) Matrix::Diagonal(x = shape) else shape
  if (is.null(scale) && is.null(shape)) {
    stop("One of scale and shape parameters must be set")
  } else if (is.null(scale)) {
    return(shape)
  } else if (is.null(shape)) {
    return(scale)
  } else {
    return(shape * scale)
  }
}

scale_and_shape_proposal <- function(sample, log_density_ratio, scale, shape) {
  scale <- scale
  shape <- shape
  list(
    sample = function(state) sample(state, get_shape_matrix(scale, shape)),
    log_density_ratio = function(state, proposed_state) log_density_ratio(
      state, proposed_state, get_shape_matrix(scale, shape)
    ),
    update = function(scale = NULL, shape = NULL) {
      if (!is.null(scale)) {
        scale <<- scale
      }
      if (!is.null(shape)) {
        shape <<- shape
      }
    },
    parameters = function() list(scale = scale, shape = shape)
  )
}
