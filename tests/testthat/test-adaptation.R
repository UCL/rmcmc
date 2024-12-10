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
    adapter, proposal, target_accept_prob, initial_scale) {
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
    adapter, dimension, check_adapter_state) {
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
        proposal <- random_walk_proposal(
          target_distribution,
          scale = 2.4 / sqrt(dimension)
        )
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
        proposal <- random_walk_proposal(target_distribution)
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
  proposal <- barker_proposal(target_distribution, scale = 1)
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
