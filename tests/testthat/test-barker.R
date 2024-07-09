barker_proposal_with_standard_normal_target <- function(step_size) {
  target_distribution <- standard_normal_target_distribution()
  barker_proposal(target_distribution, step_size = step_size)
}

for (dim in c(1, 2)) {
  test_that(
    sprintf(
      "Proposal sampling with step_size = 0 doesn't change state (dim %i)", dim
    ),
    {
      proposal <- barker_proposal_with_standard_normal_target(step_size = 0)
      withr::with_seed(seed = default_seed(), code <- {
        state <- chain_state(rnorm(dim))
        proposed_state <- proposal$sample(state)
      })
      expect_identical(proposed_state$position(), state$position())
    }
  )
  for (step_size in c(0.5, 1.)) {
    test_that(
      sprintf(
        "Proposal sampling generates valid state (dim %i, step_size %.1f)",
        dim, step_size
      ),
      {
        proposal <- barker_proposal_with_standard_normal_target(step_size)
        withr::with_seed(seed = default_seed(), code = {
          state <- chain_state(rnorm(dim))
          proposed_state <- proposal$sample(state)
        })
        check_chain_state(proposed_state)
      }
    )

    test_that(
      sprintf(
        "Proposal sampling doesn't change initial state (dim %i, step_size %.1f)",
        dim, step_size
      ),
      {
        proposal <- barker_proposal_with_standard_normal_target(step_size)
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
        "Proposal sampling changes state (dim %i, step_size %.1f)",
        dim, step_size
      ),
      {
        proposal <- barker_proposal_with_standard_normal_target(step_size)
        withr::with_seed(seed = default_seed(), code = {
          state <- chain_state(rnorm(1))
          proposed_state <- proposal$sample(state)
        })
        expect_all_different(state$position(), proposed_state$position())
      }
    )

    test_that(
      sprintf(
        "Proposal sampling with same seed gives same state (dim %i, step_size %.1f)",
        dim, step_size
      ),
      {
        proposal <- barker_proposal_with_standard_normal_target(step_size)
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
        "Proposal sampling with different seed gives different state (dim %i, step_size %.1f)",
        dim, step_size
      ),
      {
        proposal <- barker_proposal_with_standard_normal_target(step_size)
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
