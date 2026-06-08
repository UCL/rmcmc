test_that("greet works", {
  expect_output(greet("Peiyi"), regexp = "Hello Peiyi")
})
