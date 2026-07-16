test_that("as_matrix_input passes matrices through unchanged", {
  m <- matrix(1:6, nrow = 2)
  expect_identical(as_matrix_input(m), m)
})

test_that("as_matrix_input treats a bare vector as N points of 1 dimension", {
  m <- as_matrix_input(c(1, 2, 3))
  expect_equal(dim(m), c(3, 1))
  expect_equal(as.numeric(m), c(1, 2, 3))
})

test_that("as_matrix_input rejects a matrix with 0 columns", {
  expect_error(
    as_matrix_input(matrix(nrow = 3, ncol = 0)),
    "0 columns"
  )
})

test_that("map_pairs computes f for every (i, j) pair", {
  X1 <- matrix(c(0, 1, 2), ncol = 1)
  X2 <- matrix(c(0, 10), ncol = 1)
  result <- map_pairs(function(x, y) sum(x) + sum(y), X1, X2)
  expect_equal(dim(result), c(3, 2))
  expect_equal(result, matrix(c(0, 1, 2, 10, 11, 12), nrow = 3))
})

test_that("map_diag computes f only for aligned (i, i) pairs", {
  X1 <- matrix(c(0, 1, 2), ncol = 1)
  X2 <- matrix(c(10, 20, 30), ncol = 1)
  result <- map_diag(function(x, y) sum(x) + sum(y), X1, X2)
  expect_equal(result, c(10, 21, 32))
})
