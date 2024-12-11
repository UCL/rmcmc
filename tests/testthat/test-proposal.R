test_that("is_non_scalar_vector works", {
  expect_identical(is_non_scalar_vector(c(1, 2)), TRUE)
  expect_identical(is_non_scalar_vector(3.), FALSE)
  expect_identical(is_non_scalar_vector(NULL), FALSE)
  expect_identical(is_non_scalar_vector(diag(2)), FALSE)
})

test_that("get_shape_matrix works", {
  expect_identical(get_shape_matrix(1, 1), 1)
  expect_identical(get_shape_matrix(1, NULL), 1)
  expect_identical(get_shape_matrix(NULL, 1), 1)
  expect_identical(get_shape_matrix(0.5, 3), 1.5)
  expect_identical(get_shape_matrix(0.5, 3), 1.5)
  expect_identical(get_shape_matrix(2, c(3, 0.5)), c(6, 1))
  expect_identical(get_shape_matrix(NULL, c(3, 2)), c(3, 2))
  expect_identical(get_shape_matrix(0.5, diag(3, 2)), diag(1.5, 2))
  expect_identical(get_shape_matrix(NULL, diag(3)), diag(3))
  expect_error(get_shape_matrix(NULL, NULL), "must be set")
  expect_error(get_shape_matrix(-1, 2), "non-negative")
  expect_error(get_shape_matrix(c(1, 2), 1), "scalar")
})

for (dimension in c(1, 2)) {
  test_that(sprintf("scale_and_shape_proposal works in dimension %i", dimension), {
    withr::local_seed(default_seed())
    sample <- function(state, target_distribution, scale_and_shape) {
      offset <- Matrix::drop(scale_and_shape %*% rnorm(dimension))
      chain_state(position = state$position() + offset)
    }
    log_density_ratio <- function(...) 0
    proposal <- scale_and_shape_proposal(
      sample, log_density_ratio, NULL, NULL, NULL, function(d) NULL
    )
    target_distribution <- standard_normal_target_distribution()
    state <- chain_state(position = rnorm(dimension))
    expect_error(proposal$sample(state, target_distribution), "must be set")
    expect_identical(proposal$parameters(), list(scale = NULL, shape = NULL))
    proposal$update(scale = 1.)
    expect_identical(proposal$parameters(), list(scale = 1., shape = NULL))
    proposal$update(shape = 1:dimension)
    expect_identical(proposal$parameters(), list(scale = 1., shape = 1:dimension))
    proposed_state <- proposal$sample(state, target_distribution)
    check_chain_state(proposed_state)
    expect_all_different(state$position(), proposed_state$position())
    expect_identical(
      proposal$log_density(state, proposed_state, target_distribution), 0
    )
    expect_identical(proposal$default_target_accept_prob(), NULL)
    expect_identical(proposal$default_initial_scale(dimension), NULL)
  })
}
