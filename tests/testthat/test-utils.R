test_that("logistic_sigmoid works", {
  expect_equal(logistic_sigmoid(-Inf), 0)
  expect_equal(logistic_sigmoid(Inf), 1)
  expect_equal(logistic_sigmoid(0), 0.5)
  expect_equal(logistic_sigmoid(0.2), 1 - logistic_sigmoid(-0.2))
})

test_that("log1p_exp works", {
  expect_equal(log1p_exp(0.5), log1p(exp(0.5)))
  expect_equal(log1p_exp(0), log(2))
  expect_equal(log1p_exp(-Inf), 0)
})

test_that("min_1_exp works", {
  expect_equal(min_1_exp(0), 1)
  expect_equal(min_1_exp(-1.), exp(-1))
  expect_equal(min_1_exp(1), 1)
})
