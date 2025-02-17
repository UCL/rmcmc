get_shape_matrix <- function(scale, shape) {
  if (!is.null(scale) && (length(scale) > 1 || scale < 0)) {
    stop("Scale should be a non-negative scalar")
  }
  if (is.null(scale) && is.null(shape)) {
    stop("One of scale and shape parameters must be set")
  } else if (is.null(scale)) {
    shape
  } else if (is.null(shape)) {
    scale
  } else {
    shape * scale
  }
}

scale_and_shape_proposal <- function(
    sample,
    log_density_ratio,
    scale, shape,
    default_target_accept_prob,
    default_initial_scale) {
  scale <- scale
  shape <- shape
  list(
    sample = function(state, target_distribution) {
      sample(state, target_distribution, get_shape_matrix(scale, shape))
    },
    log_density_ratio = function(state, proposed_state, target_distribution) {
      shape_matrix <- get_shape_matrix(scale, shape)
      log_density_ratio(
        state, proposed_state, target_distribution, shape_matrix
      )
    },
    update = function(scale = NULL, shape = NULL) {
      if (!is.null(scale)) {
        scale <<- scale
      }
      if (!is.null(shape)) {
        shape <<- shape
      }
    },
    parameters = function() list(scale = scale, shape = shape),
    default_target_accept_prob = function() default_target_accept_prob,
    default_initial_scale = default_initial_scale
  )
}
