test_scale_and_shape_proposal(
  barker_proposal,
  proposal_name = "Barker proposal",
  target_distribution = standard_normal_target_distribution(),
  dimensions = c(1L, 2L),
  scales = c(0.5, 1., 2.)
)

test_scale_and_shape_proposal(
  bimodal_barker_proposal,
  proposal_name = "Bimodal Barker proposal",
  target_distribution = standard_normal_target_distribution(),
  dimensions = c(1L, 2L),
  scales = c(0.5, 1., 2.)
)
