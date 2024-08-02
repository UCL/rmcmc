get_scale_and_shape_from_case <- function(case, scale, dimension) {
  shape <- switch(case,
    matrix_shape_null_scale = diag(scale, dimension),
    matrix_shape_with_scale = diag(dim),
    vector_shape_null_scale = rep(scale, dimension),
    vector_shape_with_scale = rep(1, dimension),
    stop("Invalid case")
  )
  scale <- switch(case,
    matrix_shape_null_scale = ,
    vector_shape_null_scale = NULL,
    matrix_shape_with_scale = ,
    vector_shape_with_scale = scale,
    stop("Invalid case")
  )
  list(scale = scale, shape = shape)
}

test_that_no_change_when_proposal_sampling_with_zero_scale <- function(
    get_proposal_given_scale_and_shape,
    proposal_name,
    dimensions) {
  for (dimension in dimensions) {
    test_that(
      sprintf(
        "%s sampling with scale = 0 doesn't change state (dimension %i)",
        proposal_name, dimension
      ),
      {
        proposal <- get_proposal_given_scale_and_shape(scale = 0)
        withr::with_seed(seed = default_seed(), code <- {
          state <- chain_state(rnorm(dimension))
          proposed_state <- proposal$sample(state)
        })
        expect_identical(proposed_state$position(), state$position())
      }
    )
  }
}

test_that_proposal_sampling_generates_valid_state <- function(
    get_proposal_given_scale_and_shape,
    proposal_name,
    dimensions,
    scales) {
  for (dimension in dimensions) {
    for (scale in scales) {
      test_that(
        sprintf(
          "%s sampling generates valid state (dimension %i, scale %.1f)",
          proposal_name, dimension, scale
        ),
        {
          proposal <- get_proposal_given_scale_and_shape(scale)
          withr::with_seed(seed = default_seed(), code = {
            state <- chain_state(rnorm(dimension))
            proposed_state <- proposal$sample(state)
          })
          check_chain_state(proposed_state)
        }
      )
    }
  }
}

test_that_inital_state_unchanged_by_proposal_sampling <- function(
    get_proposal_given_scale_and_shape,
    proposal_name,
    dimensions,
    scales) {
  for (dimension in dimensions) {
    for (scale in scales) {
      test_that(
        sprintf(
          "%s sampling doesn't change initial state (dimension %i, scale %.1f)",
          proposal_name, dimension, scale
        ),
        {
          proposal <- get_proposal_given_scale_and_shape(scale)
          withr::with_seed(seed = default_seed(), code = {
            position <- rnorm(dimension)
            state <- chain_state(position)
            proposed_state <- proposal$sample(state)
          })
          expect_identical(state$position(), position)
        }
      )
    }
  }
}

test_that_proposal_sampling_changes_states <- function(
    get_proposal_given_scale_and_shape,
    proposal_name,
    dimensions,
    scales) {
  for (dimension in dimensions) {
    for (scale in scales) {
      test_that(
        sprintf(
          "%s sampling changes state (dimension %i, scale %.1f)",
          proposal_name, dimension, scale
        ),
        {
          proposal <- get_proposal_given_scale_and_shape(scale)
          withr::with_seed(seed = default_seed(), code = {
            state <- chain_state(rnorm(dimension))
            proposed_state <- proposal$sample(state)
          })
          expect_all_different(state$position(), proposed_state$position())
        }
      )
    }
  }
}

test_that_proposal_sampling_with_same_seed_gives_same_state <- function(
    get_proposal_given_scale_and_shape,
    proposal_name,
    dimensions,
    scales) {
  for (dimension in dimensions) {
    for (scale in scales) {
      test_that(
        sprintf(
          "%s sampling with same seed gives same state (dimension %i, scale %.1f)",
          proposal_name, dimension, scale
        ),
        {
          proposal <- get_proposal_given_scale_and_shape(scale)
          withr::with_seed(seed = default_seed(), code = {
            position <- rnorm(dimension)
            state <- chain_state(position)
            withr::with_preserve_seed({
              proposed_state <- proposal$sample(state)
            })
            proposed_state_same_seed <- proposal$sample(state)
          })
          expect_identical(
            proposed_state$position(), proposed_state_same_seed$position()
          )
        }
      )
    }
  }
}

test_that_proposal_sampling_with_different_seed_changes_state <- function(
    get_proposal_given_scale_and_shape,
    proposal_name,
    dimensions,
    scales) {
  for (dimension in dimensions) {
    for (scale in scales) {
      test_that(
        sprintf(
          "%s sampling with different seed changes state (dimension %i, scale %.1f)",
          proposal_name, dimension, scale
        ),
        {
          proposal <- get_proposal_given_scale_and_shape(scale)
          withr::with_seed(seed = default_seed(), code = {
            state <- chain_state(rnorm(dimension))
            proposed_state <- proposal$sample(state)
            proposed_state_different_seed <- proposal$sample(state)
          })
          expect_all_different(
            proposed_state$position(), proposed_state_different_seed$position()
          )
        }
      )
    }
  }
}

test_that_proposal_log_density_ratio_valid <- function(
    get_proposal_given_scale_and_shape,
    proposal_name,
    dimensions,
    scales) {
  for (dimension in dimensions) {
    for (scale in scales) {
      test_that(
        sprintf(
          "%s log density ratio valid (dimension %i, scale %.1f)",
          proposal_name, dimension, scale
        ),
        {
          proposal <- get_proposal_given_scale_and_shape(scale)
          withr::with_seed(seed = default_seed(), code = {
            state <- chain_state(rnorm(dimension))
            proposed_state <- proposal$sample(state)
          })
          log_density_ratio_forward <- proposal$log_density_ratio(state, proposed_state)
          expect_length(log_density_ratio_forward, 1)
          expect_true(is.numeric(log_density_ratio_forward))
          expect_true(is.finite(log_density_ratio_forward))
          log_density_ratio_backward <- proposal$log_density_ratio(proposed_state, state)
          expect_identical(log_density_ratio_forward, -log_density_ratio_backward)
        }
      )
    }
  }
}

test_that_proposal_with_scaled_identity_shape_equivalent_to_scale <- function(
    get_proposal_given_scale_and_shape,
    proposal_name,
    dimensions,
    scales) {
  for (dimension in dimensions) {
    for (scale in scales) {
      for (case in c("matrix_shape_null_scale", "vector_shape_null_scale")) {
        test_that(
          sprintf(
            paste0(
              "%s sampling with scaled identity shape equivalent to just ",
              "scale (dimension %i, scale %.1f, %s)"
            ),
            proposal_name, dimension, scale, case
          ),
          {
            scale_and_shape <- get_scale_and_shape_from_case(case, scale, dimension)
            proposal_scale <- get_proposal_given_scale_and_shape(scale)
            proposal_scale_and_shape <- get_proposal_given_scale_and_shape(
              scale_and_shape$scale, scale_and_shape$shape
            )
            withr::with_seed(seed = default_seed(), code = {
              position <- rnorm(dimension)
              state <- chain_state(position)
              withr::with_preserve_seed({
                proposed_state_scale <- proposal_scale$sample(state)
              })
              proposed_state_scale_and_shape <- proposal_scale_and_shape$sample(state)
            })
            expect_identical(
              proposed_state_scale$position(),
              proposed_state_scale_and_shape$position()
            )
            expect_identical(
              proposal_scale$log_density_ratio(state, proposed_state_scale),
              proposal_scale_and_shape$log_density_ratio(state, proposed_state_scale_and_shape)
            )
          }
        )
      }
    }
  }
}

test_scale_and_shape_proposal <- function(
    get_proposal_given_scale_and_shape,
    proposal_name,
    dimensions,
    scales) {
  test_that_no_change_when_proposal_sampling_with_zero_scale(
    get_proposal_given_scale_and_shape,
    proposal_name = proposal_name,
    dimensions = dimensions
  )

  test_that_proposal_sampling_generates_valid_state(
    get_proposal_given_scale_and_shape,
    proposal_name = proposal_name,
    dimensions = dimensions,
    scales = scales
  )

  test_that_inital_state_unchanged_by_proposal_sampling(
    get_proposal_given_scale_and_shape,
    proposal_name = proposal_name,
    dimensions = dimensions,
    scales = scales
  )


  test_that_proposal_sampling_changes_states(
    get_proposal_given_scale_and_shape,
    proposal_name = proposal_name,
    dimensions = dimensions,
    scales = scales
  )

  test_that_proposal_sampling_with_same_seed_gives_same_state(
    get_proposal_given_scale_and_shape,
    proposal_name = proposal_name,
    dimensions = dimensions,
    scales = scales
  )

  test_that_proposal_sampling_with_different_seed_changes_state(
    get_proposal_given_scale_and_shape,
    proposal_name = proposal_name,
    dimensions = dimensions,
    scales = scales
  )

  test_that_proposal_log_density_ratio_valid(
    get_proposal_given_scale_and_shape,
    proposal_name = proposal_name,
    dimensions = dimensions,
    scales = scales
  )

  test_that_proposal_with_scaled_identity_shape_equivalent_to_scale(
    get_proposal_given_scale_and_shape,
    proposal_name = proposal_name,
    dimensions = dimensions,
    scales = scales
  )
}
