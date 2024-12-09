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

test_that("is_non_scalar_vector works", {
  expect_true(is_non_scalar_vector(c(1, 2)))
  expect_false(is_non_scalar_vector(matrix(1:4, 2, 2)))
  expect_false(is_non_scalar_vector(3))
})

test_that("matrix vector multiply %@% op works", {
  m <- matrix(1:9, nrow = 3, ncol = 3)
  v <- c(0.5, 2, -1.2)
  s <- 0.8
  scalar_m <- matrix(2, 1, 1)
  expect_equal(drop(m %*% v), m %@% v)
  expect_equal(drop(v %*% m), v %@% m)
  expect_equal(drop(diag(v) %*% v), v %@% v)
  expect_equal(s * v, s %@% v)
  expect_equal(s * v, v %@% s)
  expect_equal(s * s, s %@% s)
  expect_equal(drop(scalar_m %*% s), scalar_m %@% s)
  expect_equal(drop(s %*% scalar_m), s %@% scalar_m)
  expect_error(m %@% s, regexp = "at least one vector")
  expect_error(s %@% m, regexp = "at least one vector")
  expect_error(m %@% m, regexp = "at least one vector")
})
