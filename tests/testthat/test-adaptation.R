check_adapter <- function(adapter) {
  expect_named(adapter, c("initialize", "update", "finalize"), ignore.order = TRUE)
  expect_type(adapter$initialize, "closure")
  expect_type(adapter$update, "closure")
  if (!is.null(adapter$finalize)) {
    expect_type(adapter$finalize, "closure")
  }
}

dummy_proposal_with_scale_parameter <- function(scale = NULL) {
  list(
    update = function(scale) scale <<- scale,
    parameters = function() list(scale = scale)
  )
}

dummy_proposal_with_shape_parameter <- function(shape = NULL) {
  list(
    update = function(shape) shape <<- shape,
    parameters = function() list(shape = shape)
  )
}

for (target_accept_prob in c(0.2, 0.4, 0.6)) {
  for (initial_scale in c(0.5, 1., 2.)) {
    for (kappa in c(0.5, 0.6, 0.8)) {
      test_that(
        sprintf(
          paste0(
            "Scale adapter works with target_accept_prob %.1f ",
            "initial_scale %.1f kappa %.1f"
          ),
          target_accept_prob, initial_scale, kappa
        ),
        {
          proposal <- dummy_proposal_with_scale_parameter()
          adapter <- scale_adapter(
            proposal = proposal,
            initial_scale = initial_scale,
            target_accept_prob = target_accept_prob,
            kappa = kappa
          )
          check_adapter(adapter)
          adapter$initialize(initial_state = NULL)
          old_scale <- initial_scale
          # If accept probability higher than target scale should be increased
          for (s in 1:2) {
            adapter$update(
              s, list(statistics = list(accept_prob = target_accept_prob + 0.1))
            )
            scale <- proposal$parameters()$scale
            expect_gt(scale, old_scale)
            old_scale <- scale
          }
          # If accept probability lower than target scale should be decreased
          for (s in 3:4) {
            adapter$update(
              s, list(statistics = list(accept_prob = target_accept_prob - 0.1))
            )
            scale <- proposal$parameters()$scale
            expect_lt(scale, old_scale)
            old_scale <- scale
          }
          # For a smooth decreasing relation between accept probability and
          # scale should adapt over long run to give close to target accept
          # probability
          for (s in 5:2000) {
            accept_prob <- exp(-scale)
            adapter$update(
              s, list(statistics = list(accept_prob = accept_prob))
            )
            scale <- proposal$parameters()$scale
          }
          expect_equal(scale, -log(target_accept_prob), tolerance = 1e-2)
        }
      )
    }
  }
}

for (dimension in c(1L, 2L, 5L)) {
  for (kappa in c(0.5, 0.6, 0.8)) {
    for (correlation in c(0.0, 0.5, 0.8)) {
      test_that(
        sprintf(
          paste0(
            "Variance adapter works for AR1 process with correlation %.1f in ",
            "dimension %i with kappa %.1f"
          ),
          correlation, dimension, kappa
        ),
        {
          proposal <- dummy_proposal_with_shape_parameter()
          adapter <- variance_adapter(proposal = proposal, kappa = kappa)
          check_adapter(adapter)
          withr::local_seed(default_seed())
          target_scales <- exp(2 * rnorm(dimension))
          position <- rnorm(dimension) * target_scales
          state <- chain_state(position = position)
          adapter$initialize(state)
          # Proposal shape parameter should be adapted to close to target scales
          # over long run
          for (s in 1:3000) {
            # Sample new position from AR1 process with stationary distribution
            # zero-mean normal with standard deviations equal to target scales
            position <- (
              correlation * position
                + sqrt(1 - correlation^2) * target_scales * rnorm(dimension)
            )
            state$update(position = position)
            adapter$update(s, list(state = state))
          }
          expect_equal(
            proposal$parameters()$shape,
            target_scales,
            tolerance = 1e-1
          )
        }
      )
    }
  }
}
