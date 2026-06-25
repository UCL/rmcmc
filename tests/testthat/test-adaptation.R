check_adapter <- function(adapter) {
  expect_named(
    adapter,
    c("initialize", "update", "finalize", "state"),
    ignore.order = TRUE
  )
  expect_type(adapter$initialize, "closure")
  expect_type(adapter$update, "closure")
  if (!is.null(adapter$finalize)) {
    expect_type(adapter$finalize, "closure")
  }
  expect_type(adapter$state, "closure")
}

dummy_proposal_with_scale_parameter <- function(scale = NULL) {
  list(
    update = function(scale) scale <<- scale,
    parameters = function() list(scale = scale),
    default_target_accept_prob = function() 0.234,
    default_initial_scale = function(dimension) 1 / sqrt(dimension)
  )
}

dummy_proposal_with_shape_parameter <- function(shape = NULL) {
  list(
    update = function(shape) shape <<- shape,
    parameters = function() list(shape = shape),
    default_target_accept_prob = function() 0.234,
    default_initial_scale = function(dimension) 1 / sqrt(dimension)
  )
}

check_scale_adapter_coerces_to_target_accept_prob <- function(
  adapter, proposal, target_accept_prob, initial_scale
) {
  # For a smooth decreasing relation between accept probability and
  # scale should adapt over long run to give close to target accept
  # probability
  scale <- initial_scale
  for (sample_index in 1:2000) {
    accept_prob <- exp(-scale)
    adapter$update(
      proposal,
      sample_index,
      list(statistics = list(accept_prob = accept_prob))
    )
    scale <- proposal$parameters()$scale
  }
  if (!is.null(adapter$finalize)) {
    adapter$finalize(proposal)
    scale <- proposal$parameters()$scale
  }
  expect_equal(scale, -log(target_accept_prob), tolerance = 1e-2)
}

check_scale_adapter_with_default_args_works <- function(
  adapter, dimension, check_adapter_state
) {
  check_adapter(adapter)
  proposal <- dummy_proposal_with_scale_parameter()
  adapter$initialize(proposal, chain_state(rep(0, dimension)))
  adapter_state <- adapter$state()
  check_adapter_state(adapter_state)
  expected_log_scale <- log(proposal$default_initial_scale(dimension))
  expect_equal(adapter_state$log_scale, expected_log_scale)
  adapter$update(
    proposal, 1, list(statistics = list(accept_prob = 1.))
  )
  adapter_state <- adapter$state()
  expect_gte(adapter_state$log_scale, expected_log_scale)
}

check_stochastic_approximation_scale_adapter_state <- function(adapter_state) {
  expect_named(adapter_state, c("log_scale"))
  expect_length(adapter_state$log_scale, 1)
}

for (target_accept_prob in c(0.2, 0.4, 0.6)) {
  for (initial_scale in c(0.5, 1., 2.)) {
    for (kappa in c(0.5, 0.6, 0.8)) {
      test_that(
        sprintf(
          paste0(
            "Stochastic approximation scale adapter works with target_accept_prob %.1f ",
            "initial_scale %.1f kappa %.1f"
          ),
          target_accept_prob, initial_scale, kappa
        ),
        {
          adapter <- scale_adapter(
            algorithm = "stochastic_approximation",
            initial_scale = initial_scale,
            target_accept_prob = target_accept_prob,
            kappa = kappa
          )
          check_adapter(adapter)
          proposal <- dummy_proposal_with_scale_parameter()
          adapter$initialize(proposal, chain_state(rep(0, dimension)))
          adapter_state <- adapter$state()
          check_stochastic_approximation_scale_adapter_state(adapter_state)
          old_scale <- initial_scale
          # If accept probability higher than target scale should be increased
          for (sample_index in 1:2) {
            adapter$update(
              proposal,
              sample_index,
              list(statistics = list(accept_prob = target_accept_prob + 0.1))
            )
            expect_type(adapter$state(), "list")
            scale <- proposal$parameters()$scale
            expect_gt(scale, old_scale)
            old_scale <- scale
          }
          # If accept probability lower than target scale should be decreased
          adapter$initialize(proposal, chain_state(rep(0, dimension)))
          old_scale <- initial_scale
          for (sample_index in 1:2) {
            adapter$update(
              proposal,
              sample_index,
              list(statistics = list(accept_prob = target_accept_prob - 0.1))
            )
            scale <- proposal$parameters()$scale
            expect_lt(scale, old_scale)
            old_scale <- scale
          }
          check_scale_adapter_coerces_to_target_accept_prob(
            adapter, proposal, target_accept_prob, initial_scale
          )
        }
      )
    }
  }
}

for (dimension in c(1L, 2L, 5L)) {
  test_that(
    sprintf(
      "Stochastic approximation scale adapter with default args works in dimension %i",
      dimension
    ),
    {
      check_scale_adapter_with_default_args_works(
        scale_adapter(algorithm = "stochastic_approximation"),
        dimension,
        check_stochastic_approximation_scale_adapter_state
      )
    }
  )
}

check_dual_averaging_scale_adapter_state <- function(adapter_state) {
  expect_named(
    adapter_state,
    c("log_scale", "smoothed_log_scale", "accept_prob_error")
  )
  expect_length(adapter_state$log_scale, 1)
  expect_length(adapter_state$smoothed_log_scale, 1)
  expect_length(adapter_state$accept_prob_error, 1)
}

for (target_accept_prob in c(0.2, 0.4, 0.6)) {
  for (initial_scale in c(0.5, 1., 2.)) {
    for (kappa in c(0.6, 0.8)) {
      for (gamma in c(0.01, 0.05)) {
        test_that(
          sprintf(
            paste0(
              "Dual averaging scale adapter works with target_accept_prob %.1f ",
              "initial_scale %.1f kappa %.1f gamma %.2f"
            ),
            target_accept_prob, initial_scale, kappa, gamma
          ),
          {
            adapter <- scale_adapter(
              algorithm = "dual_averaging",
              initial_scale = initial_scale,
              target_accept_prob = target_accept_prob,
              kappa = kappa,
              gamma = gamma
            )
            check_adapter(adapter)
            proposal <- dummy_proposal_with_scale_parameter()
            adapter$initialize(proposal, chain_state(rep(0, dimension)))
            adapter_state <- adapter$state()
            check_dual_averaging_scale_adapter_state(adapter_state)
            check_scale_adapter_coerces_to_target_accept_prob(
              adapter, proposal, target_accept_prob, initial_scale
            )
          }
        )
      }
    }
  }
}

for (dimension in c(1L, 2L, 5L)) {
  test_that(
    sprintf(
      "Dual averaging scale adapter with default args works in dimension %i",
      dimension
    ),
    {
      check_scale_adapter_with_default_args_works(
        scale_adapter(algorithm = "dual_averaging"),
        dimension,
        check_dual_averaging_scale_adapter_state
      )
    }
  )
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
          adapter <- shape_adapter(type = "variance", kappa = kappa)
          check_adapter(adapter)
          withr::local_seed(default_seed())
          target_scales <- exp(2 * rnorm(dimension))
          position <- rnorm(dimension) * target_scales
          state <- chain_state(position = position)
          adapter$initialize(proposal, state)
          adapter_state <- adapter$state()
          expect_named(
            adapter_state,
            c("mean_estimate", "variance_estimate"),
            ignore.order = TRUE
          )
          expect_length(adapter_state$mean_estimate, dimension)
          expect_length(adapter_state$variance_estimate, dimension)
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
            adapter$update(proposal, s, list(state = state))
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

for (dimension in c(1L, 2L, 3L)) {
  for (kappa in c(0.6, 1.)) {
    test_that(
      sprintf(
        "Covariance shape adapter works in dimension %i with kappa %.1f",
        dimension, kappa
      ),
      {
        withr::local_seed(default_seed())
        covariance <- random_covariance_matrix(dimension)
        chol_covariance <- t(chol(covariance))
        target_distribution <- multivariate_normal_target_distribution(
          mean = 0, covariance = covariance
        )
        proposal <- random_walk_proposal(scale = 2.4 / sqrt(dimension))
        adapter <- shape_adapter(type = "covariance", kappa = kappa)
        check_adapter(adapter)
        state <- chain_state(position = rnorm(dimension))
        adapter$initialize(proposal, state)
        adapter_state <- adapter$state()
        expect_named(
          adapter_state, c("mean_estimate", "chol_covariance_estimate")
        )
        expect_length(adapter_state$mean_estimate, dimension)
        expect_nrow(adapter_state$chol_covariance_estimate, dimension)
        expect_ncol(adapter_state$chol_covariance_estimate, dimension)
        for (sample_index in 1:10000) {
          state_and_statistics <- sample_metropolis_hastings(
            state, target_distribution, proposal
          )
          adapter$update(proposal, sample_index, state_and_statistics)
          state <- state_and_statistics$state
        }
        expect_equal(
          proposal$parameters()$shape, chol_covariance,
          tolerance = 0.1
        )
      }
    )
  }
}

for (dimension in c(1L, 2L, 3L)) {
  for (target_accept_prob in c(0.234, 0.4)) {
    test_that(
      sprintf(
        paste0(
          "Robust shape adapter works in dimension %i with ",
          "target_accept_prob %.2f "
        ),
        dimension, target_accept_prob
      ),
      {
        withr::local_seed(default_seed())
        covariance <- random_covariance_matrix(dimension)
        chol_covariance <- t(chol(covariance))
        target_distribution <- multivariate_normal_target_distribution(
          mean = 0, covariance = covariance
        )
        proposal <- random_walk_proposal()
        adapter <- robust_shape_adapter(
          kappa = 0.6,
          target_accept_prob = target_accept_prob
        )
        check_adapter(adapter)
        state <- chain_state(position = rnorm(dimension))
        adapter$initialize(proposal, state)
        adapter_state <- adapter$state()
        expect_named(adapter_state, "shape")
        expect_nrow(adapter_state$shape, dimension)
        expect_ncol(adapter_state$shape, dimension)
        mean_accept_prob <- 0.
        for (sample_index in 1:10000) {
          state_and_statistics <- sample_metropolis_hastings(
            state, target_distribution, proposal
          )
          adapter$update(proposal, sample_index, state_and_statistics)
          state <- state_and_statistics$state
          mean_accept_prob <- mean_accept_prob + (
            state_and_statistics$statistics$accept_prob - mean_accept_prob
          ) / sample_index
        }
        expect_equal(mean_accept_prob, target_accept_prob, tolerance = 0.1)
        # Proposal shape parameter should be adapted to close to Cholesky
        # factor of covariance of target distribution modulo a scaling
        # constant over long run
        in_lower_triangle <- lower.tri(chol_covariance, TRUE)
        ratios <- (
          proposal$parameters()$shape[in_lower_triangle]
          / chol_covariance[in_lower_triangle]
        )
        expect_lt(diff(range(ratios)) / median(ratios), 0.1)
      }
    )
  }
}

for (dimension in c(1L, 2L, 5L)) {
  test_that(
    sprintf(
      "Robust shape adapter default args works in dimension %i", dimension
    ),
    {
      adapter <- robust_shape_adapter()
      check_adapter(adapter)
      proposal <- dummy_proposal_with_shape_parameter()
      adapter$initialize(proposal, chain_state(rep(0, dimension)))
      adapter_state <- adapter$state()
      expect_named(adapter_state, "shape")
      initial_shape <- adapter_state$shape
      expect_nrow(initial_shape, dimension)
      expect_ncol(initial_shape, dimension)
      expect_equal(initial_shape, diag(dimension) / sqrt(dimension))
      adapter$update(
        proposal,
        1,
        list(
          proposed_state = chain_state(
            position = NULL, momentum = rep(1, dimension)
          ),
          statistics = list(accept_prob = 1.)
        )
      )
      adapter_state <- adapter$state()
      expect_gt(norm(initial_shape - adapter_state$shape), 0)
    }
  )
}

test_that("sample_chain works with dummy adapter with required interface", {
  dummy_adapter <- list(
    initialize = function(proposal, initial_state) {},
    update = function(proposal, sample_index, state_and_statistics) {},
    finalize = function(proposal) {},
    state = function() list()
  )
  target_distribution <- standard_normal_target_distribution()
  proposal <- barker_proposal(scale = 1)
  expect_no_error(
    sample_chain(
      target_distribution = target_distribution,
      proposal = proposal,
      initial_state = chain_state(0),
      n_warm_up_iteration = 1,
      n_main_iteration = 0,
      adapters = list(dummy_adapter)
    )
  )
})

test_that("Initialising scale_adapter with unrecognised algorithm gives error", {
  expect_error(scale_adapter(algorithm = "foo"), "Unrecognized algorithm")
})

test_that("Initialising shape_adapter with unrecognised type gives error", {
  expect_error(shape_adapter(type = "foo"), "Unrecognized type")
})

# ── Adapter state carry-over between stages ───────────────────────────────────

test_that(
  "stochastic_approximation_scale_adapter initialize reads current proposal scale when initial_scale not specified",
  {
    proposal <- dummy_proposal_with_scale_parameter()
    adapter <- stochastic_approximation_scale_adapter()
    # First initialise normally so the proposal acquires a default scale
    adapter$initialize(proposal, chain_state(rep(0, 2)))
    default_log_scale <- adapter$state()$log_scale
    # Manually set proposal to a different scale to simulate end of a prior stage
    carried_scale <- exp(default_log_scale) * 3.7
    proposal$update(scale = carried_scale)
    # Re-initialise (as happens at the start of a new stage): should pick up
    # the carried scale, not reset to the default
    adapter$initialize(proposal, chain_state(rep(0, 2)))
    expect_equal(adapter$state()$log_scale, log(carried_scale))
  }
)

test_that(
  "stochastic_approximation_scale_adapter explicit initial_scale overrides current proposal scale",
  {
    explicit_scale <- 5.
    proposal <- dummy_proposal_with_scale_parameter()
    adapter <- stochastic_approximation_scale_adapter(initial_scale = explicit_scale)
    # Set a different scale on the proposal to confirm it is ignored
    proposal$update(scale = 99.)
    adapter$initialize(proposal, chain_state(rep(0, 2)))
    expect_equal(adapter$state()$log_scale, log(explicit_scale))
  }
)

test_that(
  "dual_averaging_scale_adapter initialize reads current proposal scale when initial_scale not specified",
  {
    proposal <- dummy_proposal_with_scale_parameter()
    adapter <- dual_averaging_scale_adapter()
    adapter$initialize(proposal, chain_state(rep(0, 2)))
    default_log_scale <- adapter$state()$log_scale
    carried_scale <- exp(default_log_scale) * 2.5
    proposal$update(scale = carried_scale)
    adapter$initialize(proposal, chain_state(rep(0, 2)))
    expect_equal(adapter$state()$log_scale, log(carried_scale))
  }
)

test_that(
  "dual_averaging_scale_adapter explicit initial_scale overrides current proposal scale",
  {
    explicit_scale <- 3.
    proposal <- dummy_proposal_with_scale_parameter()
    adapter <- dual_averaging_scale_adapter(initial_scale = explicit_scale)
    proposal$update(scale = 99.)
    adapter$initialize(proposal, chain_state(rep(0, 2)))
    expect_equal(adapter$state()$log_scale, log(explicit_scale))
  }
)

test_that(
  "variance_shape_adapter initialize reads current proposal shape when initial_shape not specified",
  {
    proposal <- dummy_proposal_with_shape_parameter()
    adapter <- variance_shape_adapter()
    # Simulate end of a prior adaptation stage: proposal has a non-unit shape
    carried_shape <- c(2., 0.5)
    proposal$update(shape = carried_shape)
    adapter$initialize(proposal, chain_state(rep(0, 2)))
    # variance_estimate should be carried_shape^2, not rep(1, 2)
    expect_equal(adapter$state()$variance_estimate, carried_shape^2)
  }
)

test_that(
  "variance_shape_adapter explicit initial_shape overrides current proposal shape",
  {
    explicit_shape <- c(3., 0.1)
    proposal <- dummy_proposal_with_shape_parameter()
    adapter <- variance_shape_adapter(initial_shape = explicit_shape)
    # Set a different shape on the proposal to confirm it is ignored
    proposal$update(shape = c(99., 99.))
    adapter$initialize(proposal, chain_state(rep(0, 2)))
    expect_equal(adapter$state()$variance_estimate, explicit_shape^2)
  }
)

test_that(
  "variance_shape_adapter falls back to unit variances when no current proposal shape and no initial_shape",
  {
    proposal <- dummy_proposal_with_shape_parameter() # shape starts as NULL
    adapter <- variance_shape_adapter()
    adapter$initialize(proposal, chain_state(rep(0, 2)))
    expect_equal(adapter$state()$variance_estimate, rep(1., 2))
  }
)

test_that(
  "variance_shape_adapter falls back to unit variances when current proposal shape is matrix-valued",
  {
    # A matrix-valued shape (from covariance adapter) cannot be used as a
    # diagonal variance initialisation — should fall back to identity
    proposal <- dummy_proposal_with_shape_parameter()
    proposal$update(shape = diag(2)) # matrix, not a vector
    adapter <- variance_shape_adapter()
    adapter$initialize(proposal, chain_state(rep(0, 2)))
    expect_equal(adapter$state()$variance_estimate, rep(1., 2))
  }
)

test_that(
  "covariance_shape_adapter initialize reads current proposal shape when initial_shape not specified",
  {
    proposal <- dummy_proposal_with_shape_parameter()
    adapter <- covariance_shape_adapter()
    # Simulate end of a prior stage: proposal has a non-identity Cholesky shape
    carried_shape <- matrix(c(2., 0., 0.5, 1.), nrow = 2)
    proposal$update(shape = carried_shape)
    adapter$initialize(proposal, chain_state(rep(0, 2)))
    expect_equal(adapter$state()$chol_covariance_estimate, carried_shape)
  }
)

test_that(
  "covariance_shape_adapter explicit initial_shape overrides current proposal shape",
  {
    explicit_shape <- matrix(c(3., 0., 0.2, 0.8), nrow = 2)
    proposal <- dummy_proposal_with_shape_parameter()
    adapter <- covariance_shape_adapter(initial_shape = explicit_shape)
    proposal$update(shape = diag(99., 2))
    adapter$initialize(proposal, chain_state(rep(0, 2)))
    expect_equal(adapter$state()$chol_covariance_estimate, explicit_shape)
  }
)

test_that(
  "covariance_shape_adapter falls back to identity when no current proposal shape and no initial_shape",
  {
    proposal <- dummy_proposal_with_shape_parameter() # shape starts as NULL
    adapter <- covariance_shape_adapter()
    adapter$initialize(proposal, chain_state(rep(0, 2)))
    expect_equal(adapter$state()$chol_covariance_estimate, diag(1., 2))
  }
)

test_that(
  "covariance_shape_adapter falls back to identity when current proposal shape is a vector",
  {
    # A vector shape (from variance adapter) cannot be used as a Cholesky
    # factor — should fall back to identity
    proposal <- dummy_proposal_with_shape_parameter()
    proposal$update(shape = c(2., 0.5)) # vector, not a matrix
    adapter <- covariance_shape_adapter()
    adapter$initialize(proposal, chain_state(rep(0, 2)))
    expect_equal(adapter$state()$chol_covariance_estimate, diag(1., 2))
  }
)

test_that(
  "scale adapter log_scale does not reset to default at stage boundary in two-stage sample_chain",
  {
    # Run a two-stage chain where the first stage has enough iterations for the
    # scale adapter to move away from its default. Check that warm_up_statistics
    # shows the log_scale does NOT jump back to the default value at the stage
    # boundary (iteration 51 vs 50).
    target_distribution <- standard_normal_target_distribution()
    dimension <- 2L
    default_log_scale <- log(
      barker_proposal()$default_initial_scale(dimension)
    )
    s1 <- stochastic_approximation_scale_adapter(initial_scale = 1.)
    s2 <- stochastic_approximation_scale_adapter() # no explicit initial_scale
    withr::with_seed(default_seed(), {
      results <- sample_chain(
        target_distribution = target_distribution,
        initial_state = rnorm(dimension),
        n_warm_up_iteration = 100L,
        n_main_iteration = 1L,
        adapters = list(list(s1, 50L), list(s2)),
        trace_warm_up = TRUE,
        show_progress_bar = FALSE
      )
    })
    log_scale_col <- results$warm_up_statistics[, "log_scale"]
    # At the stage boundary (row 51) the log_scale should be close to row 50,
    # not jump back to the default value
    jump_at_boundary <- abs(log_scale_col[51] - log_scale_col[50])
    reset_to_default <- abs(log_scale_col[51] - default_log_scale)
    # If carry-over works, jump_at_boundary << reset_to_default
    expect_lt(jump_at_boundary, reset_to_default)
  }
)

# ── progressive_adaptation_schedule tests ─────────────────────────────────────

test_that("progressive_adaptation_schedule returns a function", {
  schedule <- progressive_adaptation_schedule()
  expect_type(schedule, "closure")
})

test_that(
  "progressive_adaptation_schedule with n_warm_up_iteration > n_fixed + n_diagonal produces three stages",
  {
    # Default n_fixed=50, n_diagonal=50; 200 > 100 so dense stage gets remainder
    schedule <- progressive_adaptation_schedule(
      n_fixed_shape_iteration = 50L,
      n_diagonal_shape_iteration = 50L
    )
    stages <- schedule(200L)
    expect_length(stages, 3L)
  }
)

test_that(
  "progressive_adaptation_schedule stage iteration counts are correct when three stages produced",
  {
    schedule <- progressive_adaptation_schedule(
      n_fixed_shape_iteration = 50L,
      n_diagonal_shape_iteration = 50L
    )
    stages <- schedule(200L)
    # Extract the trailing integer from each stage spec
    iter_counts <- sapply(stages, function(s) s[[length(s)]])
    expect_equal(iter_counts[[1]], 50L)
    expect_equal(iter_counts[[2]], 50L)
    # Last stage has no trailing integer (it is the remainder stage passed
    # as a no-count stage to check_and_process_adapters), so verify via
    # parsing the full staged list
    parsed <- rmcmc:::check_and_process_adapters(stages, 200L)
    expect_equal(parsed[[1]]$n_iteration, 50L)
    expect_equal(parsed[[2]]$n_iteration, 50L)
    expect_equal(parsed[[3]]$n_iteration, 100L)
  }
)

test_that(
  paste0(
    "progressive_adaptation_schedule with n_warm_up_iteration = n_fixed + n_diagonal ",
    "produces two stages with no dense stage"
  ),
  {
    # n_dense = 100 - 50 - 50 = 0 -> dense stage is skipped
    schedule <- progressive_adaptation_schedule(
      n_fixed_shape_iteration = 50L,
      n_diagonal_shape_iteration = 50L
    )
    stages <- schedule(100L)
    expect_length(stages, 2L)
    parsed <- rmcmc:::check_and_process_adapters(stages, 100L)
    expect_equal(parsed[[1]]$n_iteration, 50L)
    expect_equal(parsed[[2]]$n_iteration, 50L)
  }
)

test_that(
  paste0(
    "progressive_adaptation_schedule with n_warm_up_iteration < n_fixed + n_diagonal ",
    "falls back to two stages splitting iterations equally"
  ),
  {
    # 60 < 50 + 50 = 100, so fallback: stage1 = 60 %/% 2 = 30, stage2 = 30
    schedule <- progressive_adaptation_schedule(
      n_fixed_shape_iteration = 50L,
      n_diagonal_shape_iteration = 50L
    )
    stages <- schedule(60L)
    expect_length(stages, 2L)
    parsed <- rmcmc:::check_and_process_adapters(stages, 60L)
    expect_equal(parsed[[1]]$n_iteration, 30L)
    expect_equal(parsed[[2]]$n_iteration, 30L)
  }
)

test_that(
  "progressive_adaptation_schedule with n_warm_up_iteration = 1 does not error",
  {
    schedule <- progressive_adaptation_schedule()
    expect_no_error(schedule(1L))
  }
)

test_that(
  "progressive_adaptation_schedule respects custom n_fixed_shape_iteration and n_diagonal_shape_iteration",
  {
    schedule <- progressive_adaptation_schedule(
      n_fixed_shape_iteration = 20L,
      n_diagonal_shape_iteration = 30L
    )
    parsed <- rmcmc:::check_and_process_adapters(schedule(100L), 100L)
    expect_equal(parsed[[1]]$n_iteration, 20L)
    expect_equal(parsed[[2]]$n_iteration, 30L)
    expect_equal(parsed[[3]]$n_iteration, 50L)
  }
)

test_that(
  "progressive_adaptation_schedule stage 1 contains only a scale adapter",
  {
    schedule <- progressive_adaptation_schedule()
    stages <- schedule(200L)
    parsed <- rmcmc:::check_and_process_adapters(stages, 200L)
    # Stage 1 should have exactly one adapter (the scale adapter)
    expect_length(parsed[[1]]$adapters, 1L)
    expect_true(rmcmc:::is_adapter(parsed[[1]]$adapters[[1]]))
  }
)

test_that(
  "progressive_adaptation_schedule stage 2 contains scale and diagonal shape adapters",
  {
    schedule <- progressive_adaptation_schedule()
    stages <- schedule(200L)
    parsed <- rmcmc:::check_and_process_adapters(stages, 200L)
    # Stage 2 should have exactly two adapters
    expect_length(parsed[[2]]$adapters, 2L)
    expect_true(rmcmc:::is_adapter(parsed[[2]]$adapters[[1]]))
    expect_true(rmcmc:::is_adapter(parsed[[2]]$adapters[[2]]))
  }
)

test_that(
  "progressive_adaptation_schedule stage 3 contains scale and dense shape adapters",
  {
    schedule <- progressive_adaptation_schedule()
    stages <- schedule(200L)
    parsed <- rmcmc:::check_and_process_adapters(stages, 200L)
    # Stage 3 should have exactly two adapters
    expect_length(parsed[[3]]$adapters, 2L)
    expect_true(rmcmc:::is_adapter(parsed[[3]]$adapters[[1]]))
    expect_true(rmcmc:::is_adapter(parsed[[3]]$adapters[[2]]))
  }
)

test_that(
  "progressive_adaptation_schedule result passes through check_and_process_adapters without error",
  {
    for (n in c(1L, 50L, 100L, 200L, 500L)) {
      schedule <- progressive_adaptation_schedule()
      expect_no_error(
        rmcmc:::check_and_process_adapters(schedule(n), n)
      )
    }
  }
)

test_that(
  "sample_chain with progressive_adaptation_schedule() runs without error and returns correct structure",
  {
    target_distribution <- standard_normal_target_distribution()
    n_warm_up <- 120L
    n_main <- 10L
    withr::with_seed(default_seed(), {
      results <- sample_chain(
        target_distribution = target_distribution,
        initial_state = rnorm(2),
        n_warm_up_iteration = n_warm_up,
        n_main_iteration = n_main,
        adapters = progressive_adaptation_schedule(
          n_fixed_shape_iteration = 40L,
          n_diagonal_shape_iteration = 40L
        ),
        trace_warm_up = FALSE,
        show_progress_bar = FALSE
      )
    })
    expect_named(
      results,
      c("final_state", "traces", "statistics"),
      ignore.order = TRUE
    )
    expect_nrow(results$traces, n_main)
    expect_nrow(results$statistics, n_main)
  }
)
