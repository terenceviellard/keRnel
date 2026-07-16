# General-purpose numeric linear algebra helpers for the (quasi-)symmetric
# positive-definite matrices evaluate() produces -- not tied to any
# particular downstream model. Motivated by a recurring pattern across
# consumers of keRnel covariance matrices (see dev/ameliorations/README.md,
# item 1d): base R's chol2inv() drops dimnames, and a near-singular matrix
# makes chol() fail outright with no recovery path.

#' Robust Cholesky factor (adaptive jitter), dimnames preserved
#'
#' @description
#' Like [chol()], but with a recovery path: a (quasi-)singular `mat` no
#' longer makes it fail outright -- a small jitter is added to the
#' diagonal, retried with a 10x larger jitter each time `chol()` still
#' fails, until it succeeds. `dimnames()` are preserved either way.
#'
#' @param mat A symmetric (quasi-)positive-definite numeric matrix.
#' @param pen_diag Initial jitter added to the diagonal; multiplied by 10 on
#'   each `chol()` failure.
#' @return The upper-triangular factor `R` such that `mat + jitter = t(R) %*% R`,
#'   with the same `dimnames()` as `mat` and an attribute `pen_diag_used`
#'   (the jitter that was actually retained).
#' @export
#' @examples
#' k <- se_kernel(length_scale = 1)
#' K <- evaluate(k, c(0, 1, 2))
#' robust_chol(K)
robust_chol <- function(mat, pen_diag = 1e-10) {
  nm <- dimnames(mat)
  diag(mat) <- diag(mat) + pen_diag
  R <- tryCatch(
    chol(mat),
    error = function(e) {
      diag(mat) <- diag(mat) - pen_diag
      robust_chol(mat, 10 * pen_diag)
    }
  )
  if (is.null(attr(R, "pen_diag_used"))) attr(R, "pen_diag_used") <- pen_diag
  dimnames(R) <- nm
  R
}

#' Robust matrix inverse via Cholesky, dimnames preserved
#'
#' @inheritParams robust_chol
#' @return The inverse of `mat + jitter` (see [robust_chol()]), with the
#'   same `dimnames()` as `mat`.
#' @export
#' @examples
#' k <- se_kernel(length_scale = 1)
#' K <- evaluate(k, c(0, 1, 2))
#' robust_chol_inv(K)
robust_chol_inv <- function(mat, pen_diag = 1e-10) {
  nm <- dimnames(mat)
  inv <- chol2inv(robust_chol(mat, pen_diag))
  dimnames(inv) <- nm
  inv
}
