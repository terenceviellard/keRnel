# These tests are the safety net for the performance work in
# R/distance-matrix.R: every vectorised fast path must agree exactly (up to
# floating-point tolerance) with the always-correct, slower reference
# (map_pairs() applied to the scalar distance/pairwise function).

test_that("squared_euclidean distance_matrix agrees with the scalar formula", {
  X1 <- matrix(c(0, 1, 2, 0, 1, 4), ncol = 2)
  X2 <- matrix(c(0, 3, 0, -1), ncol = 2)
  expect_equal(
    distance_matrix(squared_euclidean_distance, X1, X2),
    map_pairs(squared_euclidean_distance, X1, X2)
  )
})

test_that("euclidean distance_matrix agrees with the scalar formula", {
  X1 <- matrix(c(0, 1, 2, 0, 1, 4), ncol = 2)
  X2 <- matrix(c(0, 3, 0, -1), ncol = 2)
  expect_equal(
    distance_matrix(euclidean_distance, X1, X2),
    map_pairs(euclidean_distance, X1, X2)
  )
})

test_that("manhattan distance_matrix agrees with the scalar formula", {
  X1 <- matrix(c(0, 1, 2, 0, 1, 4), ncol = 2)
  X2 <- matrix(c(0, 3, 0, -1), ncol = 2)
  expect_equal(
    distance_matrix(manhattan_distance, X1, X2),
    map_pairs(manhattan_distance, X1, X2)
  )
})

test_that("dot_product distance_matrix agrees with the scalar formula", {
  X1 <- matrix(c(0, 1, 2, 0, 1, 4), ncol = 2)
  X2 <- matrix(c(0, 3, 0, -1), ncol = 2)
  expect_equal(
    distance_matrix(dot_product, X1, X2),
    map_pairs(dot_product, X1, X2)
  )
})

test_that("equality distance_matrix agrees with the scalar formula", {
  X1 <- matrix(c(0, 1, 2, 0, 1, 4), ncol = 2)
  X2 <- matrix(c(0, 3, 0, -1, 2, 4), ncol = 2)
  expect_equal(
    distance_matrix(equality, X1, X2),
    map_pairs(equality, X1, X2)
  )
})

test_that("distance_matrix falls back to map_pairs() for an unrecognised distance function", {
  X1 <- matrix(c(0, 1, 2), ncol = 1)
  X2 <- matrix(c(0, 3), ncol = 1)
  expect_equal(
    distance_matrix(chebyshev_distance, X1, X2),
    map_pairs(chebyshev_distance, X1, X2)
  )
})

test_that("squared_euclidean distance_matrix never returns a negative value", {
  X1 <- matrix(rep(1, 4), ncol = 1)
  X2 <- matrix(rep(1, 4), ncol = 1)
  expect_true(all(distance_matrix(squared_euclidean_distance, X1, X2) >= 0))
})

test_that("pairwise_matrix.stationary_kernel agrees exactly with the scalar pairwise() reference", {
  k <- se_kernel(length_scale = 1.3)
  X1 <- matrix(c(-1, 0, 0.5, 2), ncol = 1)
  X2 <- matrix(c(0, 1, 3), ncol = 1)

  fast <- pairwise_matrix(k, X1, X2)
  reference <- map_pairs(function(x, y) pairwise(k, x, y), X1, X2)

  expect_equal(fast, reference, tolerance = 1e-10)
})

test_that("pairwise_matrix.stationary_kernel also agrees on multi-dimensional inputs", {
  k <- se_kernel(length_scale = 0.7)
  X1 <- matrix(c(0, 0, 1, 2, -1, 3), ncol = 2, byrow = TRUE)
  X2 <- matrix(c(0, 0, 2, 2), ncol = 2, byrow = TRUE)

  fast <- pairwise_matrix(k, X1, X2)
  reference <- map_pairs(function(x, y) pairwise(k, x, y), X1, X2)

  expect_equal(fast, reference, tolerance = 1e-10)
})

test_that("evaluate() on a stationary kernel uses the fast path transparently", {
  k <- se_kernel(length_scale = 1)
  x <- c(-2, -1, 0, 1, 2)
  expect_equal(
    evaluate(k, x),
    map_pairs(function(a, b) pairwise(k, a, b), as_matrix_input(x), as_matrix_input(x)),
    tolerance = 1e-10
  )
})

# Regression test: a custom stationary_kernel/dot_product_kernel that only
# defines pairwise() (no shape()) is a perfectly valid, if slower, kernel --
# it must still work through evaluate()'s default dense_engine path rather
# than crashing with "no applicable method for 'shape'" (pairwise_matrix's
# fast path used to call shape() unconditionally). See
# vignette("02-intermediate")'s "Writing your own kernel" section, which
# relies on exactly this working.

test_that("a stationary_kernel without a shape() method still evaluates correctly", {
  cosine_kernel <- function(length_scale) {
    structure(
      list(
        length_scale = wrapped_param(length_scale),
        distance_function = euclidean_distance,
        engine = dense_engine()
      ),
      class = c("cosine_kernel_no_shape", "stationary_kernel", "kernel")
    )
  }
  pairwise.cosine_kernel_no_shape <- function(kernel, x1, x2) {
    cos(kernel$distance_function(x1, x2) / unwrap_param(kernel$length_scale))
  }
  # Register the S3 method locally for this test (normally done via @export
  # in a package, see R/distance-matrix.R's has_s3_method()).
  registerS3method("pairwise", "cosine_kernel_no_shape", pairwise.cosine_kernel_no_shape, envir = environment())

  k <- cosine_kernel(length_scale = 1)
  x <- c(0, 1, 2)
  expect_equal(
    evaluate(k, x),
    map_pairs(function(a, b) pairwise(k, a, b), as_matrix_input(x), as_matrix_input(x)),
    tolerance = 1e-10
  )
})

test_that("a dot_product_kernel without a shape() method still evaluates correctly", {
  triple_kernel <- function(slope_var) {
    structure(
      list(slope_var = wrapped_param(slope_var), distance_function = dot_product, engine = dense_engine()),
      class = c("triple_kernel_no_shape", "dot_product_kernel", "kernel")
    )
  }
  pairwise.triple_kernel_no_shape <- function(kernel, x1, x2) {
    unwrap_param(kernel$slope_var) * kernel$distance_function(x1, x2)^3
  }
  registerS3method("pairwise", "triple_kernel_no_shape", pairwise.triple_kernel_no_shape, envir = environment())

  k <- triple_kernel(slope_var = 1)
  x <- c(0, 1, 2)
  expect_equal(
    evaluate(k, x),
    map_pairs(function(a, b) pairwise(k, a, b), as_matrix_input(x), as_matrix_input(x)),
    tolerance = 1e-10
  )
})

test_that("has_s3_method() correctly detects an existing method", {
  expect_true(has_s3_method("shape", se_kernel(length_scale = 1)))
  expect_true(has_s3_method("shape", linear_kernel(slope_var = 1)))
})

test_that("has_s3_method() correctly detects a missing method", {
  broken <- structure(list(), class = c("no_such_class_anywhere", "kernel"))
  expect_false(has_s3_method("shape", broken))
})

# diag_distance()/pairwise_diag(): the diagonal-only counterparts backing
# diagonal_engine() (R/engines.R). Same safety net as above, one dimension
# lower.

test_that("squared_euclidean diag_distance agrees with the scalar formula", {
  X1 <- matrix(c(0, 1, 2, 0, 1, 4), ncol = 2)
  X2 <- matrix(c(0, 3, 1, 0, -1, 2), ncol = 2)
  expect_equal(
    diag_distance(squared_euclidean_distance, X1, X2),
    map_diag(squared_euclidean_distance, X1, X2)
  )
})

test_that("euclidean diag_distance agrees with the scalar formula", {
  X1 <- matrix(c(0, 1, 2, 0, 1, 4), ncol = 2)
  X2 <- matrix(c(0, 3, 1, 0, -1, 2), ncol = 2)
  expect_equal(
    diag_distance(euclidean_distance, X1, X2),
    map_diag(euclidean_distance, X1, X2)
  )
})

test_that("manhattan diag_distance agrees with the scalar formula", {
  X1 <- matrix(c(0, 1, 2, 0, 1, 4), ncol = 2)
  X2 <- matrix(c(0, 3, 1, 0, -1, 2), ncol = 2)
  expect_equal(
    diag_distance(manhattan_distance, X1, X2),
    map_diag(manhattan_distance, X1, X2)
  )
})

test_that("dot_product diag_distance agrees with the scalar formula", {
  X1 <- matrix(c(0, 1, 2, 0, 1, 4), ncol = 2)
  X2 <- matrix(c(0, 3, 1, 0, -1, 2), ncol = 2)
  expect_equal(
    diag_distance(dot_product, X1, X2),
    map_diag(dot_product, X1, X2)
  )
})

test_that("equality diag_distance agrees with the scalar formula", {
  X1 <- matrix(c(0, 1, 2, 0, 1, 4), ncol = 2)
  X2 <- matrix(c(0, 3, 2, 0, 1, 4), ncol = 2)
  expect_equal(
    diag_distance(equality, X1, X2),
    map_diag(equality, X1, X2)
  )
})

test_that("diag_distance falls back to map_diag() for an unrecognised distance function", {
  X1 <- matrix(c(0, 1, 2), ncol = 1)
  X2 <- matrix(c(0, 3, -1), ncol = 1)
  expect_equal(
    diag_distance(chebyshev_distance, X1, X2),
    map_diag(chebyshev_distance, X1, X2)
  )
})

test_that("pairwise_diag.stationary_kernel agrees exactly with the scalar pairwise() reference", {
  k <- se_kernel(length_scale = 1.3)
  X1 <- matrix(c(-1, 0, 0.5, 2), ncol = 1)
  X2 <- matrix(c(0, 1, -2, 3), ncol = 1)

  fast <- pairwise_diag(k, X1, X2)
  reference <- map_diag(function(x, y) pairwise(k, x, y), X1, X2)

  expect_equal(fast, reference, tolerance = 1e-10)
})

test_that("pairwise_diag.dot_product_kernel agrees exactly with the scalar pairwise() reference", {
  k <- linear_kernel(slope_var = 1.7)
  X1 <- matrix(c(0, 0, 1, 2, -1, 3), ncol = 2, byrow = TRUE)
  X2 <- matrix(c(0, 1, 2, -1, 1, 0), ncol = 2, byrow = TRUE)

  fast <- pairwise_diag(k, X1, X2)
  reference <- map_diag(function(x, y) pairwise(k, x, y), X1, X2)

  expect_equal(fast, reference, tolerance = 1e-10)
})
