# constant_kernel: k(x1, x2) = value, for any inputs

test_that("constant_kernel returns the same value regardless of inputs", {
  k <- constant_kernel(value = 3)
  expect_equal(pairwise(k, c(1, 2), c(3, 4)), 3)
  expect_equal(pairwise(k, c(0, 0), c(0, 0)), 3)
})

test_that("constant_kernel accepts negative values", {
  k <- constant_kernel(value = -2)
  expect_equal(pairwise(k, c(1), c(2)), -2)
})

test_that("constant_kernel defaults to a value of 1", {
  k <- constant_kernel()
  expect_equal(pairwise(k, c(1), c(2)), 1)
})

test_that("constant_kernel's evaluate() fills the whole matrix with value", {
  k <- constant_kernel(value = 4)
  K <- evaluate(k, c(1, 2, 3))
  expect_equal(K, matrix(4, nrow = 3, ncol = 3))
})

test_that("constant_kernel's value is optimisable (wrapped_param)", {
  k <- constant_kernel(value = 4)
  expect_equal(unname(get_params(k)), 4, tolerance = 1e-10)
})

test_that("constant_kernel rejects non-finite/non-numeric values", {
  expect_error(constant_kernel(value = NA), "finite numeric")
  expect_error(constant_kernel(value = NaN), "finite numeric")
  expect_error(constant_kernel(value = Inf), "finite numeric")
  expect_error(constant_kernel(value = "x"), "finite numeric")
  expect_error(constant_kernel(value = NULL), "finite numeric")
})

# variance_kernel: k(x1, x2) = variance, for any inputs, variance > 0

test_that("variance_kernel returns the same positive value regardless of inputs", {
  k <- variance_kernel(variance = 2.5)
  expect_equal(pairwise(k, c(1, 2), c(3, 4)), 2.5)
})

test_that("variance_kernel defaults to a variance of 1", {
  k <- variance_kernel()
  expect_equal(pairwise(k, c(1), c(2)), 1)
})

test_that("variance_kernel's evaluate() fills the whole matrix with variance", {
  k <- variance_kernel(variance = 7)
  K <- evaluate(k, c(1, 2, 3))
  expect_equal(K, matrix(7, nrow = 3, ncol = 3))
})

test_that("variance_kernel rejects non-positive variance", {
  expect_error(variance_kernel(variance = 0), "strictly positive")
  expect_error(variance_kernel(variance = -1), "strictly positive")
})

test_that("variance_kernel's variance is optimisable", {
  k <- variance_kernel(variance = 2)
  expect_equal(unname(get_params(k)), 2, tolerance = 1e-10)
})
