barker_proposal_with_standard_normal_target <- function(
    scale = NULL, shape = NULL) {
  target_distribution <- standard_normal_target_distribution()
  barker_proposal(target_distribution, scale = scale, shape = shape)
}

for (dim in c(1, 2)) {
  test_that(
    sprintf(
      "Proposal sampling with scale = 0 doesn't change state (dim %i)", dim
    ),
    {
      proposal <- barker_proposal_with_standard_normal_target(scale = 0)
      withr::with_seed(seed = default_seed(), code <- {
        state <- chain_state(rnorm(dim))
        proposed_state <- proposal$sample(state)
      })
      expect_identical(proposed_state$position(), state$position())
    }
  )
  for (scale in c(0.5, 1.)) {
    test_that(
      sprintf(
        "Proposal sampling generates valid state (dim %i, scale %.1f)",
        dim, scale
      ),
      {
        proposal <- barker_proposal_with_standard_normal_target(scale)
        withr::with_seed(seed = default_seed(), code = {
          state <- chain_state(rnorm(dim))
          proposed_state <- proposal$sample(state)
        })
        check_chain_state(proposed_state)
      }
    )

    test_that(
      sprintf(
        "Proposal sampling doesn't change initial state (dim %i, scale %.1f)",
        dim, scale
      ),
      {
        proposal <- barker_proposal_with_standard_normal_target(scale)
        withr::with_seed(seed = default_seed(), code = {
          position <- rnorm(dim)
          state <- chain_state(position)
          proposed_state <- proposal$sample(state)
        })
        expect_identical(state$position(), position)
      }
    )

    test_that(
      sprintf(
        "Proposal sampling changes state (dim %i, scale %.1f)",
        dim, scale
      ),
      {
        proposal <- barker_proposal_with_standard_normal_target(scale)
        withr::with_seed(seed = default_seed(), code = {
          state <- chain_state(rnorm(1))
          proposed_state <- proposal$sample(state)
        })
        expect_all_different(state$position(), proposed_state$position())
      }
    )

    test_that(
      sprintf(
        "Proposal sampling with same seed gives same state (dim %i, scale %.1f)",
        dim, scale
      ),
      {
        proposal <- barker_proposal_with_standard_normal_target(scale)
        withr::with_seed(seed = default_seed(), code = {
          position <- rnorm(dim)
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

    test_that(
      sprintf(
        "Proposal sampling with different seed gives different state (dim %i, scale %.1f)",
        dim, scale
      ),
      {
        proposal <- barker_proposal_with_standard_normal_target(scale)
        withr::with_seed(seed = default_seed(), code = {
          state <- chain_state(rnorm(dim))
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


for (dim in c(1, 2)) {
  for (scale in c(1., 2., 0.5)) {
    for (case in c("matrix_shape_null_scale", "vector_shape_null_scale")) {
      test_that(
        sprintf(
          paste0(
            "Proposal sampling with scaled identity shape equivalent to just ",
            "scale (dim %i, scale %.1f, %s)"
          ),
          dim, scale, case
        ),
        {
          shape <- switch(case,
            matrix_shape_null_scale = diag(scale, dim),
            matrix_shape_with_scale = diag(dim),
            vector_shape_null_scale = rep(scale, dim),
            vector_shape_with_scale = rep(1, dim),
            stop("Invalid case")
          )
          scale_for_use_with_shape <- switch(case,
            matrix_shape_null_scale = ,
            vector_shape_null_scale = NULL,
            matrix_shape_with_scale = ,
            vector_shape_with_scale = scale,
            stop("Invalid case")
          )
          proposal_scale <- barker_proposal_with_standard_normal_target(scale)
          proposal_scale_and_shape <- barker_proposal_with_standard_normal_target(
            scale_for_use_with_shape, shape
          )
          withr::with_seed(seed = default_seed(), code = {
            position <- rnorm(dim)
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
