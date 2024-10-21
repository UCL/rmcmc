for (n_warm_up_iteration in c(0, 1, 10)) {
  for (n_main_iteration in c(0, 1, 10)) {
    for (dimension in c(1, 2)) {
      for (trace_warm_up in c(TRUE, FALSE)) {
        for (show_progress_bar in c(TRUE, FALSE)) {
          for (wrapped_initial_state in c(TRUE, FALSE)) {
            test_that(
              sprintf(
                paste0(
                  "Sampling chain with %i warm-up iterations, ",
                  "%i main iterations, dimension %i, ",
                  "trace_warm_up = %i, show_progress_bar = %i ",
                  "and wrapped_initial_state = %i works"
                ),
                n_warm_up_iteration,
                n_main_iteration,
                dimension,
                trace_warm_up,
                show_progress_bar,
                wrapped_initial_state
              ),
              {
                target_distribution <- standard_normal_target_distribution()
                barker_proposal(target_distribution)
                proposal <- barker_proposal(target_distribution)
                adapters <- list(
                  scale_adapter(proposal, initial_scale = 1.)
                )
                withr::with_seed(default_seed(), {
                  position <- rnorm(dimension)
                })
                if (wrapped_initial_state) {
                  initial_state <- chain_state(position)
                } else {
                  initial_state <- position
                }
                results <- sample_chain(
                  target_distribution = target_distribution,
                  proposal = proposal,
                  initial_state = initial_state,
                  n_warm_up_iteration = n_warm_up_iteration,
                  n_main_iteration = n_main_iteration,
                  adapters = adapters,
                  trace_warm_up = trace_warm_up,
                  show_progress_bar = show_progress_bar
                )
                expected_results_names <- c(
                  "final_state", "traces", "statistics"
                )
                if (trace_warm_up) {
                  expected_results_names <- c(
                    expected_results_names,
                    "warm_up_traces",
                    "warm_up_statistics"
                  )
                }
                expect_named(
                  results,
                  expected_results_names,
                  ignore.order = TRUE,
                )
                expect_nrow(results$traces, n_main_iteration)
                expect_ncol(results$traces, dimension + 1)
                expect_nrow(results$statistics, n_main_iteration)
                expect_ncol(results$statistics, 1)
                if (trace_warm_up) {
                  expect_nrow(results$warm_up_traces, n_warm_up_iteration)
                  expect_ncol(results$warm_up_traces, dimension + 1)
                  expect_nrow(results$warm_up_statistics, n_warm_up_iteration)
                  expect_ncol(results$warm_up_statistics, 2)
                }
              }
            )
          }
        }
      }
    }
  }
}
