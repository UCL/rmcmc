test_scale_and_shape_proposal(
  langevin_proposal,
  proposal_name = "Langevin proposal",
  target_distribution = standard_normal_target_distribution(),
  dimensions = c(1L, 2L),
  scales = c(0.5, 1., 2.)
)

for (dimension in c(1L, 2L, 3L)) {
  for (scale in c(0.5, 1., 2.)) {
    test_that(
      sprintf(
        "Langevin involution is an involution (dimension %i, scale %.1f)",
        dimension, scale
      ),
      {
        target_distribution <- standard_normal_target_distribution()
        withr::with_seed(seed = default_seed(), {
          state <- chain_state(
            position = rnorm(dimension), momentum = rnorm(dimension)
          )
        })
        inv_state <- involution_langevin(state, scale, target_distribution)
        inv_inv_state <- involution_langevin(inv_state, scale, target_distribution)
        expect_equal(state$position(), inv_inv_state$position())
        expect_equal(state$momentum(), inv_inv_state$momentum())
      }
    )
  }
}
