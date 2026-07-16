#' Computation engines
#'
#' @description
#' An engine picks which mapping primitive (layer 2, see [map_pairs()]) is
#' used to turn a kernel's [pairwise()] formula into a full covariance
#' matrix. This is the *policy* layer: adding a new engine means choosing/
#' combining existing mapping primitives, not writing a new loop.
#'
#' `dense_engine()` is the default for every kernel, computing the full
#' N1 x N2 matrix via [map_pairs()].
#'
#' @param engine An engine object.
#' @param kernel A kernel object.
#' @param X1,X2 Numeric matrices (see [as_matrix_input()]).
#'
#' @return A numeric matrix (or, for other engines, a vector/array).
#' @name engine
NULL

#' @rdname engine
#' @export
engine_eval <- function(engine, kernel, X1, X2) {
  UseMethod("engine_eval")
}

#' Dense computation engine
#'
#' @description
#' Computes the full N1 x N2 covariance matrix via [map_pairs()]. The
#' default engine for every kernel.
#'
#' `dense_engine()` has no compilation constraints to work around -- it is
#' simply the natural default.
#'
#' @return An object of class `dense_engine`/`engine`.
#' @export
#' @examples
#' k <- se_kernel(length_scale = 1)
#' evaluate(k, c(1, 2, 3))
dense_engine <- function() {
  structure(list(), class = c("dense_engine", "engine"))
}

#' @export
engine_eval.dense_engine <- function(engine, kernel, X1, X2) {
  pairwise_matrix(kernel, X1, X2)
}

#' Diagonal computation engine
#'
#' @description
#' Computes only the diagonal of the covariance matrix, i.e. `k(x1[i], x2[i])`
#' for each aligned pair of points, as a vector of length `nrow(X1)` --
#' `O(n)` instead of `O(n^2)`. Useful when only the diagonal is needed, e.g.
#' the predictive variance of a Gaussian Process.
#'
#' `diagonal_engine()` is a deliberate addition: the optimisation it
#' captures (avoiding an `O(n^2)` computation when only the diagonal is
#' needed) is purely algebraic, independent of any particular compilation
#' or hardware backend.
#'
#' Requires `x1` and `x2` to have the same number of points.
#'
#' @return An object of class `diagonal_engine`/`engine`.
#' @export
#' @examples
#' k <- se_kernel(length_scale = 1, engine = diagonal_engine())
#' evaluate(k, c(0, 1, 2), c(0, 1, 2))
diagonal_engine <- function() {
  structure(list(), class = c("diagonal_engine", "engine"))
}

#' @export
engine_eval.diagonal_engine <- function(engine, kernel, X1, X2) {
  if (nrow(X1) != nrow(X2)) {
    stop(
      "diagonal_engine() requires `x1` and `x2` to have the same number of points ",
      "(got ", nrow(X1), " and ", nrow(X2), ").",
      call. = FALSE
    )
  }
  pairwise_diag(kernel, X1, X2)
}
