# Reference values below are derived by hand from each kernel's formula
# (see dev/ameliorations/design/03-analyse-risques.Rmd, risk #5: tests must be independently
# verifiable).

# linear_kernel: k(x1, x2) = slope_var * dot(x1, x2)

test_that("linear_kernel matches the closed-form value for a known dot product", {
  k <- linear_kernel(slope_var = 2)
  # dot(c(1,2), c(3,4)) == 11
  expect_equal(pairwise(k, c(1, 2), c(3, 4)), 2 * 11)
})

test_that("linear_kernel passes through the origin", {
  k <- linear_kernel(slope_var = 5)
  expect_equal(pairwise(k, c(0, 0), c(3, 4)), 0)
})

test_that("linear_kernel's evaluate() agrees with the vectorised fast path", {
  k <- linear_kernel(slope_var = 1.7)
  X <- matrix(c(0, 0, 1, 2, -1, 3), ncol = 2, byrow = TRUE)
  expect_equal(
    evaluate(k, X),
    map_pairs(function(a, b) pairwise(k, a, b), X, X),
    tolerance = 1e-10
  )
})

test_that("linear_kernel rejects non-positive slope_var", {
  expect_error(linear_kernel(slope_var = 0), "strictly positive")
})

# affine_kernel: k(x1, x2) = slope_var * dot(x1 - offset, x2 - offset)

test_that("affine_kernel matches the closed-form value for a known dot product", {
  k <- affine_kernel(slope_var = 2, offset = c(1, 1))
  # (c(1,2)-c(1,1)) . (c(3,4)-c(1,1)) = c(0,1) . c(2,3) = 3
  expect_equal(pairwise(k, c(1, 2), c(3, 4)), 2 * 3)
})

test_that("affine_kernel passes through (offset, 0)", {
  k <- affine_kernel(slope_var = 5, offset = c(2, -1))
  expect_equal(pairwise(k, c(2, -1), c(3, 4)), 0)
})

test_that("affine_kernel's evaluate() agrees with the vectorised fast path", {
  k <- affine_kernel(slope_var = 0.9, offset = c(0.5, -0.5))
  X <- matrix(c(0, 0, 1, 2, -1, 3), ncol = 2, byrow = TRUE)
  expect_equal(
    evaluate(k, X),
    map_pairs(function(a, b) pairwise(k, a, b), X, X),
    tolerance = 1e-10
  )
})

# polynomial_kernel: k(x1, x2) = (gamma * dot(x1, x2) + constant)^degree

test_that("polynomial_kernel matches the closed-form value for a known dot product", {
  k <- polynomial_kernel(degree = 2, gamma = 1, constant = 1)
  # dot(c(1,2), c(3,4)) == 11; k = (1*11 + 1)^2 = 144
  expect_equal(pairwise(k, c(1, 2), c(3, 4)), 144)
})

test_that("polynomial_kernel's evaluate() agrees with the vectorised fast path", {
  k <- polynomial_kernel(degree = 3, gamma = 0.5, constant = -1)
  X <- matrix(c(0, 0, 1, 2, -1, 3), ncol = 2, byrow = TRUE)
  expect_equal(
    evaluate(k, X),
    map_pairs(function(a, b) pairwise(k, a, b), X, X),
    tolerance = 1e-10
  )
})

test_that("polynomial_kernel rejects a non-positive-integer degree", {
  expect_error(polynomial_kernel(degree = 0, gamma = 1, constant = 0), "positive integer")
  expect_error(polynomial_kernel(degree = 1.5, gamma = 1, constant = 0), "positive integer")
})

test_that("polynomial_kernel rejects non-positive gamma", {
  expect_error(polynomial_kernel(degree = 2, gamma = 0, constant = 0), "strictly positive")
})

test_that("polynomial_kernel's degree is not mutable via kupdate()", {
  k <- polynomial_kernel(degree = 2, gamma = 1, constant = 0)
  expect_error(kupdate(k, degree = 3), "No such hyperparameter")
})

# sigmoid_kernel: k(x1, x2) = tanh(alpha * dot(x1, x2) + constant)

test_that("sigmoid_kernel matches the closed-form value for a known dot product", {
  k <- sigmoid_kernel(alpha = 1, constant = 0)
  # dot(c(1,0), c(1,0)) == 1; k = tanh(1)
  expect_equal(pairwise(k, c(1, 0), c(1, 0)), tanh(1), tolerance = 1e-10)
})

test_that("sigmoid_kernel's evaluate() agrees with the vectorised fast path", {
  k <- sigmoid_kernel(alpha = 0.3, constant = 0.2)
  X <- matrix(c(0, 0, 1, 2, -1, 3), ncol = 2, byrow = TRUE)
  expect_equal(
    evaluate(k, X),
    map_pairs(function(a, b) pairwise(k, a, b), X, X),
    tolerance = 1e-10
  )
})

test_that("sigmoid_kernel is bounded within [-1, 1]", {
  k <- sigmoid_kernel(alpha = 5, constant = 0)
  K <- evaluate(k, c(-10, -1, 0, 1, 10))
  expect_true(all(K >= -1 & K <= 1))
})

test_that("dot-product kernels carry the dot_product_kernel class", {
  expect_true(inherits(linear_kernel(slope_var = 1), "dot_product_kernel"))
  expect_true(inherits(affine_kernel(slope_var = 1, offset = 0), "dot_product_kernel"))
  expect_true(inherits(polynomial_kernel(degree = 1, gamma = 1, constant = 0), "dot_product_kernel"))
  expect_true(inherits(sigmoid_kernel(alpha = 1, constant = 0), "dot_product_kernel"))
})
