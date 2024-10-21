for (dimension in c(1L, 2L)) {
  test_that(sprintf("Creating chain state with only position in dimension %i works", dimension), {
    position <- rep(0, dimension)
    state <- chain_state(position)
    check_chain_state(state)
    expect_identical(state$position(), position)
    expect_identical(state$momentum(), NULL)
    expect_identical(state$dimension(), dimension)
  })
  test_that(sprintf("Creating chain state with momentum in dimension %i works", dimension), {
    withr::with_seed(default_seed(), {
      position <- rnorm(dimension)
      momentum <- rnorm(dimension)
    })
    state <- chain_state(position, momentum)
    check_chain_state(state)
    expect_identical(state$position(), position)
    expect_identical(state$momentum(), momentum)
    expect_identical(state$dimension(), dimension)
  })
  for (use_value_and_gradient in c(TRUE, FALSE)) {
    test_that(
      sprintf(
        paste0(
          "Evaluating log density (gradient) with chain state and ",
          "use_value_and_gradient = %i in dimension %i works"
        ),
        use_value_and_gradient,
        dimension
      ),
      {
        withr::with_seed(default_seed(), position <- rnorm(dimension))
        target_distribution <- standard_normal_target_distribution(
          use_value_and_gradient
        )
        state <- chain_state(position)
        expect_identical(
          state$log_density(target_distribution),
          target_distribution$log_density(position)
        )
        expect_identical(
          state$gradient_log_density(target_distribution),
          target_distribution$gradient_log_density(position)
        )
      }
    )
  }
  test_that(sprintf("Copying chain state in dimension %i works", dimension), {
    withr::with_seed(default_seed(), {
      position <- rnorm(dimension)
      momentum <- rnorm(dimension)
    })
    state <- chain_state(position, momentum)
    copied_state <- state$copy()
    check_chain_state(copied_state)
    expect_identical(state$position(), position)
    expect_identical(state$position(), copied_state$position())
    expect_identical(state$momentum(), momentum)
    expect_identical(state$momentum(), copied_state$momentum())
    expect_identical(state$dimension(), dimension)
    expect_identical(state$dimension(), copied_state$dimension())
  })
  test_that(sprintf("Updating chain state in dimension %i works", dimension), {
    withr::with_seed(default_seed(), {
      position <- rnorm(dimension)
      momentum <- rnorm(dimension)
      new_position <- rnorm(dimension)
      new_momentum <- rnorm(dimension)
    })
    state <- chain_state(position, momentum)
    original_state <- state$copy()
    state$update(position = new_position)
    expect_identical(original_state$position(), position)
    expect_identical(state$position(), new_position)
    state$update(momentum = new_momentum)
    expect_identical(original_state$momentum(), momentum)
    expect_identical(state$momentum(), new_momentum)
  })
  test_that(
    sprintf(
      "Evaluating log density (gradient) on updated chain state in dimension %i works",
      dimension
    ),
    {
      withr::with_seed(default_seed(), {
        position <- rnorm(dimension)
        new_position <- rnorm(dimension)
        new_position_2 <- rnorm(dimension)
      })
      state <- chain_state(position)
      target_distribution <- standard_normal_target_distribution()
      expect_identical(
        state$log_density(target_distribution),
        target_distribution$log_density(position)
      )
      expect_identical(
        state$gradient_log_density(target_distribution),
        target_distribution$gradient_log_density(position)
      )
      state$update(position = new_position)
      expect_identical(
        state$log_density(target_distribution),
        target_distribution$log_density(new_position)
      )
      expect_identical(
        state$gradient_log_density(target_distribution),
        target_distribution$gradient_log_density(new_position)
      )
      copied_state <- state$copy()
      copied_state$update(position = new_position_2)
      expect_identical(
        copied_state$log_density(target_distribution),
        target_distribution$log_density(new_position_2)
      )
      expect_identical(
        copied_state$gradient_log_density(target_distribution),
        target_distribution$gradient_log_density(new_position_2)
      )
      expect_identical(
        state$log_density(target_distribution),
        target_distribution$log_density(new_position)
      )
      expect_identical(
        state$gradient_log_density(target_distribution),
        target_distribution$gradient_log_density(new_position)
      )
    }
  )
}
