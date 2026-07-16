# kernel_interactions(): an ergonomic helper for building the additive +
# multiplicative interaction terms of a set of kernels (see
# R/kernel-interactions.R).

test_that("kernel_interactions() of a single kernel is just that kernel", {
  k1 <- se_kernel(length_scale = 1)
  k <- kernel_interactions(list(k1))
  expect_equal(evaluate(k, c(0, 1)), evaluate(k1, c(0, 1)))
})

test_that("kernel_interactions() of two kernels gives k1 + k2 + k1*k2", {
  k1 <- se_kernel(length_scale = 1)
  k2 <- linear_kernel(slope_var = 1)
  k <- kernel_interactions(list(k1, k2))

  x <- c(0, 1, 2)
  expected <- evaluate(k1, x) + evaluate(k2, x) + evaluate(k1, x) * evaluate(k2, x)
  expect_equal(evaluate(k, x), expected, tolerance = 1e-10)
})

test_that("kernel_interactions() of three kernels gives the full 7-term expansion", {
  k1 <- se_kernel(length_scale = 1)
  k2 <- linear_kernel(slope_var = 1)
  k3 <- periodic_kernel(length_scale = 1, period = 2)
  k <- kernel_interactions(list(k1, k2, k3))

  x <- c(0, 1, 2)
  e1 <- evaluate(k1, x)
  e2 <- evaluate(k2, x)
  e3 <- evaluate(k3, x)
  expected <- e1 + e2 + e3 + e1 * e2 + e1 * e3 + e2 * e3 + e1 * e2 * e3
  expect_equal(evaluate(k, x), expected, tolerance = 1e-10)
})

test_that("max_order limits the expansion to lower-order interactions", {
  k1 <- se_kernel(length_scale = 1)
  k2 <- linear_kernel(slope_var = 1)
  k3 <- periodic_kernel(length_scale = 1, period = 2)
  k <- kernel_interactions(list(k1, k2, k3), max_order = 2)

  x <- c(0, 1, 2)
  e1 <- evaluate(k1, x)
  e2 <- evaluate(k2, x)
  e3 <- evaluate(k3, x)
  # No triple product k1*k2*k3 included at max_order = 2.
  expected <- e1 + e2 + e3 + e1 * e2 + e1 * e3 + e2 * e3
  expect_equal(evaluate(k, x), expected, tolerance = 1e-10)
})

test_that("kernel_interactions() rejects an empty kernel list", {
  expect_error(kernel_interactions(list()), "at least one")
})

test_that("kernel_interactions() rejects an out-of-range max_order", {
  k1 <- se_kernel(length_scale = 1)
  expect_error(kernel_interactions(list(k1), max_order = 0), "between 1 and")
  expect_error(kernel_interactions(list(k1), max_order = 2), "between 1 and")
})
