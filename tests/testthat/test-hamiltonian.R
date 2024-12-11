hamiltonian_proposal_with_n_step <- function(n_step, sample_n_step = NULL) {
  function(scale = NULL, shape = NULL) {
    hamiltonian_proposal(
      n_step = n_step,
      scale = scale,
      shape = shape,
      sample_n_step = sample_n_step
    )
  }
}

for (n_step in c(1L, 2L, 5L)) {
  test_scale_and_shape_proposal(
    hamiltonian_proposal_with_n_step(n_step),
    proposal_name = sprintf("Hamiltonian proposal with n_step = %i", n_step),
    target_distribution = standard_normal_target_distribution(),
    dimensions = c(1L, 2L),
    scales = c(0.5, 1., 2.)
  )
}

for (dimension in c(1L, 2L, 3L)) {
  for (scale in c(0.5, 1., 2.)) {
    for (n_step in c(1L, 2L, 5L)) {
      test_that(
        sprintf(
          paste0(
            "Hamiltonian involution is an involution ",
            "(dimension %i, scale %.1f, n_step %i)"
          ),
          dimension, scale, n_step
        ),
        {
          target_distribution <- standard_normal_target_distribution()
          withr::with_seed(seed = default_seed(), {
            state <- chain_state(
              position = rnorm(dimension), momentum = rnorm(dimension)
            )
          })
          inv_state <- involution_hamiltonian(
            state, n_step, scale, target_distribution
          )
          inv_inv_state <- involution_hamiltonian(
            inv_state, n_step, scale, target_distribution
          )
          expect_equal(state$position(), inv_inv_state$position())
          expect_equal(state$momentum(), inv_inv_state$momentum())
        }
      )
    }
  }
}

test_scale_and_shape_proposal(
  hamiltonian_proposal_with_n_step(
    n_step = c(1, 5),
    sample_n_step = function(n) n[1] + sample.int(n[2] - n[1] + 1, 1) - 1
  ),
  proposal_name = sprintf("Hamiltonian proposal with randomized n_step"),
  target_distribution = standard_normal_target_distribution(),
  dimensions = c(1L, 2L),
  scales = c(0.5, 1., 2.)
)
