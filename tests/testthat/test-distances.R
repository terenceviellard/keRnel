test_that("euclidean and squared euclidean distances are consistent (3-4-5 triangle)", {
  x1 <- c(0, 0)
  x2 <- c(3, 4)
  expect_equal(euclidean_distance(x1, x2), 5)
  expect_equal(squared_euclidean_distance(x1, x2), 25)
})

test_that("manhattan and chebyshev distances on a simple example", {
  x1 <- c(0, 0)
  x2 <- c(3, 4)
  expect_equal(manhattan_distance(x1, x2), 7)
  expect_equal(chebyshev_distance(x1, x2), 4)
})

test_that("cosine distance is zero for colinear vectors and one for orthogonal vectors", {
  expect_equal(cosine_distance(c(1, 1), c(2, 2)), 0, tolerance = 1e-8)
  expect_equal(cosine_distance(c(1, 0), c(0, 1)), 1, tolerance = 1e-8)
})

test_that("minkowski distance reduces to manhattan (p=1) and euclidean (p=2)", {
  x1 <- c(0, 0)
  x2 <- c(3, 4)
  expect_equal(minkowski_distance(x1, x2, p = 1), manhattan_distance(x1, x2))
  expect_equal(minkowski_distance(x1, x2, p = 2), euclidean_distance(x1, x2))
})

test_that("hamming distance counts mismatches", {
  expect_equal(hamming_distance(c(1, 0, 1), c(1, 1, 1)), 1)
  expect_equal(hamming_distance(c(1, 1, 1), c(1, 1, 1)), 0)
})

test_that("jaccard distance is zero for identical binary vectors", {
  expect_equal(jaccard_distance(c(1, 0, 1), c(1, 0, 1)), 0)
})

test_that("mahalanobis distance reduces to euclidean distance for an identity covariance", {
  x1 <- c(0, 0)
  x2 <- c(3, 4)
  expect_equal(mahalanobis_distance(x1, x2, diag(2)), euclidean_distance(x1, x2))
})

test_that("dot_product computes the inner product", {
  expect_equal(dot_product(c(1, 2, 3), c(4, 5, 6)), 32)
})

test_that("equality returns 1 for identical vectors and 0 otherwise", {
  expect_equal(equality(c(1, 2), c(1, 2)), 1)
  expect_equal(equality(c(1, 2), c(1, 3)), 0)
})

test_that("canberra, bray_curtis and correlation distances behave on simple cases", {
  expect_equal(canberra_distance(c(1, 1), c(1, 1)), 0, tolerance = 1e-6)
  expect_equal(bray_curtis_distance(c(1, 1), c(1, 1)), 0, tolerance = 1e-6)
  expect_equal(correlation_distance(c(1, 2, 3), c(2, 4, 6)), 0, tolerance = 1e-6)
})

test_that("resolve_distance_function accepts functions and shortcut names", {
  expect_identical(resolve_distance_function(euclidean_distance), euclidean_distance)
  expect_identical(resolve_distance_function("manhattan"), manhattan_distance)
})

test_that("resolve_distance_function rejects unknown shortcuts and non-function values", {
  expect_error(resolve_distance_function("not-a-distance"), "must be a function")
  expect_error(resolve_distance_function(42), "must be a function")
})
