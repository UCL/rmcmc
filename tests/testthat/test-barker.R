barker_proposal_with_standard_normal_target <- function(
    scale = NULL, shape = NULL) {
  target_distribution <- standard_normal_target_distribution()
  barker_proposal(target_distribution, scale = scale, shape = shape)
}

test_scale_and_shape_proposal(
  barker_proposal_with_standard_normal_target,
  proposal_name = "Barker proposal",
  dimensions = c(1L, 2L),
  scales = c(0.5, 1., 2.)
)
