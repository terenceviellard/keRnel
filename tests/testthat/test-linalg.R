test_that("robust_chol() matches base chol() for a well-conditioned matrix", {
  k <- se_kernel(length_scale = 1)
  K <- evaluate(k, c(0, 1, 2))
  expect_equal(robust_chol(K, pen_diag = 0), chol(K), tolerance = 1e-8, ignore_attr = TRUE)
})

test_that("robust_chol() preserves the input's dimnames", {
  K <- matrix(c(2, 0, 0, 2), nrow = 2, dimnames = list(c("a", "b"), c("a", "b")))
  R <- robust_chol(K)
  expect_equal(dimnames(R), dimnames(K))
})

test_that("robust_chol() recovers from a singular matrix via adaptive jitter", {
  K <- matrix(c(1, 1, 1, 1), nrow = 2) # exactly singular
  R <- robust_chol(K, pen_diag = 1e-10)
  expect_true(all(is.finite(R)))
  expect_gt(attr(R, "pen_diag_used"), 0)
})

test_that("robust_chol_inv() returns the correct inverse with dimnames preserved", {
  K <- matrix(c(2, 0, 0, 2), nrow = 2, dimnames = list(c("a", "b"), c("a", "b")))
  inv <- robust_chol_inv(K, pen_diag = 0)
  expect_equal(inv %*% K, diag(2), ignore_attr = TRUE, tolerance = 1e-8)
  expect_equal(dimnames(inv), dimnames(K))
})
