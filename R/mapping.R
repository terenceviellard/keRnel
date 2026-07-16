# Layer 2 of the architecture (dev/ameliorations/design/01-architecture.Rmd): pure vectorisation
# mechanics, agnostic of both the kernel and the engine. This is the single
# place where any future performance work (better vectorisation, Rcpp...)
# should happen -- no engine or kernel should ever need to change because
# map_pairs() got faster.
#
# The current implementation is a naive, obviously-correct double loop.
# Performance is deliberately deferred (see dev/ameliorations/design/03-analyse-risques.Rmd,
# risk #4): get the maths right first, then measure, then optimise only if
# needed.

#' Coerce a kernel input to a matrix of points
#'
#' @description
#' `keRnel` always represents a collection of input points as an N x D
#' matrix (N points, D dimensions/features each). A bare numeric vector of
#' length N is treated as a convenience shorthand for N points of a single
#' dimension (`matrix(x, ncol = 1)`).
#'
#' This is a deliberate simplification: a single multi-dimensional point
#' must be passed as a 1-row matrix, avoiding any ambiguity between "N
#' points of 1 dimension" and "a single point of N dimensions".
#'
#' @param x A numeric vector or matrix.
#' @return A numeric matrix.
#' @keywords internal
as_matrix_input <- function(x) {
  if (is.matrix(x)) {
    if (ncol(x) == 0) {
      stop(
        "Input has 0 columns -- a matrix of points must have at least one ",
        "feature/dimension.",
        call. = FALSE
      )
    }
    return(x)
  }
  matrix(x, ncol = 1)
}

#' Apply a pairwise function to every pair of rows of two matrices
#'
#' @description
#' The "all pairs" mapping scheme: returns an N1 x N2 matrix whose `(i, j)`
#' entry is `f(X1[i, ], X2[j, ])`. Used by [dense_engine()].
#'
#' @param f A function of two numeric vectors (one row of `X1`, one row of
#'   `X2`), returning a single numeric value.
#' @param X1,X2 Numeric matrices, with the same number of columns.
#' @return A numeric matrix of dimension `nrow(X1)` x `nrow(X2)`.
#' @keywords internal
map_pairs <- function(f, X1, X2) {
  n1 <- nrow(X1)
  n2 <- nrow(X2)
  result <- matrix(NA_real_, nrow = n1, ncol = n2)
  for (i in seq_len(n1)) {
    xi <- X1[i, ]
    for (j in seq_len(n2)) {
      result[i, j] <- f(xi, X2[j, ])
    }
  }
  result
}

#' Apply a pairwise function to aligned rows of two matrices
#'
#' @description
#' The "aligned pairs" mapping scheme: returns a vector of length N whose
#' `i`-th entry is `f(X1[i, ], X2[i, ])`. Used by [diagonal_engine()].
#'
#' @param f A function of two numeric vectors, returning a single numeric
#'   value.
#' @param X1,X2 Numeric matrices with the same number of rows and columns.
#' @return A numeric vector of length `nrow(X1)`.
#' @keywords internal
map_diag <- function(f, X1, X2) {
  n <- nrow(X1)
  vapply(seq_len(n), function(i) f(X1[i, ], X2[i, ]), numeric(1))
}
