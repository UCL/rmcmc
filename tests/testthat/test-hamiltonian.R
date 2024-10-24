hamiltonian_proposal_with_standard_normal_target <- function(n_step) {
  function(scale = NULL, shape = NULL) {
    hamiltonian_proposal(
      target_distribution = standard_normal_target_distribution(),
      n_step = n_step,
      scale = scale,
      shape = shape
    )
  }
}

for (n_step in c(1L, 2L, 5L)) {
  test_scale_and_shape_proposal(
    hamiltonian_proposal_with_standard_normal_target(n_step),
    proposal_name = sprintf("Hamiltonian proposal with n_step = %i", n_step),
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
