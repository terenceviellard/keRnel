test_that("dense_engine computes the same result as a manual map_pairs call", {
  k <- se_kernel(length_scale = 1)
  X1 <- matrix(c(0, 1, 2), ncol = 1)
  X2 <- matrix(c(0, 1), ncol = 1)

  via_engine <- engine_eval(k$engine, k, X1, X2)
  via_map_pairs <- map_pairs(function(x, y) pairwise(k, x, y), X1, X2)

  expect_equal(via_engine, via_map_pairs)
})

test_that("engine_eval dispatches on the class of the engine, not the kernel", {
  expect_true(inherits(dense_engine(), c("dense_engine", "engine")))
})

# diagonal_engine: a specialised engine (see R/engines.R) that computes only
# the diagonal, O(n) instead of O(n^2).

test_that("diagonal_engine computes the same result as a manual map_diag call", {
  k <- se_kernel(length_scale = 1)
  X1 <- matrix(c(0, 1, 2), ncol = 1)
  X2 <- matrix(c(0, 1, 3), ncol = 1)

  via_engine <- engine_eval(diagonal_engine(), k, X1, X2)
  via_map_diag <- map_diag(function(x, y) pairwise(k, x, y), X1, X2)

  expect_equal(via_engine, via_map_diag, tolerance = 1e-10)
})

test_that("diagonal_engine matches the diagonal of dense_engine's full matrix", {
  k <- se_kernel(length_scale = 1.4)
  X <- matrix(c(0, 1, 2, -1), ncol = 1)

  full <- engine_eval(dense_engine(), k, X, X)
  diag_only <- engine_eval(diagonal_engine(), k, X, X)

  expect_equal(diag_only, diag(full), tolerance = 1e-10)
})

test_that("diagonal_engine errors when x1 and x2 have different lengths", {
  k <- se_kernel(length_scale = 1)
  X1 <- matrix(c(0, 1, 2), ncol = 1)
  X2 <- matrix(c(0, 1), ncol = 1)
  expect_error(engine_eval(diagonal_engine(), k, X1, X2), "same number of points")
})

test_that("evaluate() with diagonal_engine returns a vector, not a matrix", {
  k <- se_kernel(length_scale = 1, engine = diagonal_engine())
  result <- evaluate(k, c(0, 1, 2), c(0, 1, 3))
  expect_false(is.matrix(result))
  expect_length(result, 3)
})

test_that("evaluate() with diagonal_engine propagates NA to aligned pairs", {
  k <- se_kernel(length_scale = 1, engine = diagonal_engine())
  result <- evaluate(k, c(0, NA, 2), c(0, 1, NA))
  expect_false(is.na(result[1]))
  expect_true(is.na(result[2]))
  expect_true(is.na(result[3]))
})
