# Vectorised distance matrices, the key performance lever identified after
# measuring map_pairs() at realistic sizes (see dev/ameliorations/design/03-analyse-risques.Rmd,
# risk #4 -- the naive double loop took ~23s for a 1000x1000 matrix).
#
# Insight: for a stationary kernel, pairwise(kernel, x1, x2) is always
# `shape(kernel, distance_function(x1, x2))`, where `shape()` is an
# elementwise scalar transform (exp(), ^, etc.) that R already vectorises
# for free over a matrix. The only genuinely expensive part is computing the
# N1 x N2 distance matrix itself -- and for the handful of distances most
# kernels actually use, that has a closed-form vectorised expression that
# avoids an explicit R-level loop over points entirely.
#
# Distances without a known vectorised form fall back to map_pairs(), which
# remains correct (if slow) for any distance function, including
# user-supplied custom ones.

#' Compute a full distance matrix, vectorised when possible
#'
#' @param distance_function A distance function (see [distances]).
#' @param X1,X2 Numeric matrices.
#' @return A numeric matrix of dimension `nrow(X1)` x `nrow(X2)`.
#' @keywords internal
distance_matrix <- function(distance_function, X1, X2) {
  if (identical(distance_function, squared_euclidean_distance)) {
    return(.squared_euclidean_distance_matrix(X1, X2))
  }
  if (identical(distance_function, euclidean_distance)) {
    return(sqrt(.squared_euclidean_distance_matrix(X1, X2)))
  }
  if (identical(distance_function, manhattan_distance)) {
    return(.manhattan_distance_matrix(X1, X2))
  }
  if (identical(distance_function, dot_product)) {
    return(X1 %*% t(X2))
  }
  if (identical(distance_function, equality)) {
    return(.equality_matrix(X1, X2))
  }
  # No known vectorised form: fall back to the generic, slower mapping.
  map_pairs(distance_function, X1, X2)
}

.equality_matrix <- function(X1, X2) {
  result <- matrix(1, nrow = nrow(X1), ncol = nrow(X2))
  for (k in seq_len(ncol(X1))) {
    result <- result * outer(X1[, k], X2[, k], "==")
  }
  result
}

.squared_euclidean_distance_matrix <- function(X1, X2) {
  sq1 <- rowSums(X1^2)
  sq2 <- rowSums(X2^2)
  d2 <- outer(sq1, sq2, "+") - 2 * (X1 %*% t(X2))
  # Clamp tiny negative values caused by floating-point cancellation
  # (mathematically d2 >= 0 always).
  pmax(d2, 0)
}

.manhattan_distance_matrix <- function(X1, X2) {
  d <- matrix(0, nrow = nrow(X1), ncol = nrow(X2))
  for (k in seq_len(ncol(X1))) {
    d <- d + abs(outer(X1[, k], X2[, k], "-"))
  }
  d
}

#' Pairwise covariance matrix, vectorised when the kernel supports it
#'
#' @description
#' Default implementation: falls back to the always-correct, slower
#' [map_pairs()] applied to [pairwise()]. Stationary kernels get a much
#' faster method (`pairwise_matrix.stationary_kernel`) via
#' [distance_matrix()] + [shape()].
#'
#' @param kernel A kernel object.
#' @param X1,X2 Numeric matrices.
#' @return A numeric matrix of dimension `nrow(X1)` x `nrow(X2)`.
#' @keywords internal
pairwise_matrix <- function(kernel, X1, X2) {
  UseMethod("pairwise_matrix")
}

#' @keywords internal
#' @exportS3Method
pairwise_matrix.default <- function(kernel, X1, X2) {
  map_pairs(function(x, y) pairwise(kernel, x, y), X1, X2)
}

# A stationary_kernel/dot_product_kernel only gets the fast shape()-based
# path if it actually provides a shape() method -- a custom kernel that
# defines pairwise() but not shape() (a perfectly valid, if slower, kernel:
# see vignette("b-intermediate")'s "Writing your own kernel" section) must
# still work correctly through the default dense_engine, just without the
# speed-up. Checked once per call via has_s3_method() rather than assumed,
# so a missing shape() falls back to the always-correct map_pairs()/
# map_diag() path instead of erroring with "no applicable method".

#' @keywords internal
has_s3_method <- function(generic, x) {
  any(vapply(class(x), function(cls) !is.null(utils::getS3method(generic, cls, optional = TRUE)), logical(1)))
}

#' @keywords internal
#' @exportS3Method
pairwise_matrix.stationary_kernel <- function(kernel, X1, X2) {
  if (!has_s3_method("shape", kernel)) {
    return(pairwise_matrix.default(kernel, X1, X2))
  }
  d <- distance_matrix(kernel$distance_function, X1, X2)
  shape(kernel, d)
}

#' @keywords internal
#' @exportS3Method
pairwise_matrix.dot_product_kernel <- function(kernel, X1, X2) {
  if (!has_s3_method("shape", kernel)) {
    return(pairwise_matrix.default(kernel, X1, X2))
  }
  dp <- distance_matrix(kernel$distance_function, X1, X2)
  shape(kernel, dp)
}

#' Elementwise shape function of a stationary kernel
#'
#' @description
#' For a stationary kernel, `pairwise(kernel, x1, x2)` is always
#' `shape(kernel, distance_function(x1, x2))`: a scalar function of the
#' distance alone. Expressing it this way lets [pairwise_matrix()] apply it
#' to a whole distance matrix at once -- R vectorises elementwise math
#' (`exp()`, `^`, ...) over a matrix for free, so this is the one piece each
#' stationary kernel must provide to get a fast, vectorised
#' [evaluate()] for free, on top of its scalar [pairwise()] (which remains
#' the reference definition of correctness, exercised directly in tests).
#'
#' A custom kernel's `shape()` method is entirely optional -- a
#' `stationary_kernel`/`dot_product_kernel` that only defines [pairwise()]
#' remains correct, just without this speed-up (see
#' `vignette("b-intermediate")`, "Writing your own kernel", and
#' `vignette("c-expert")`, "Vectorisation internals").
#'
#' @param kernel A kernel object.
#' @param d A numeric value, vector, or matrix of distances.
#' @return The same shape as `d`.
#' @export
#' @examples
#' k <- se_kernel(length_scale = 1)
#' shape(k, c(0, 1, 4)) # exp(-0.5 * d / length_scale^2)
shape <- function(kernel, d) {
  UseMethod("shape")
}

# The diagonal-only counterpart of distance_matrix()/pairwise_matrix(),
# backing diagonal_engine() (R/engines.R). Same idea, one dimension lower:
# an N-length vector of aligned-pair distances/dot-products instead of an
# N1 x N2 matrix, computed via rowSums() rather than outer().

#' Compute a vector of point-aligned distances, vectorised when possible
#'
#' @param distance_function A distance function (see [distances]).
#' @param X1,X2 Numeric matrices with the same number of rows.
#' @return A numeric vector of length `nrow(X1)`.
#' @keywords internal
diag_distance <- function(distance_function, X1, X2) {
  if (identical(distance_function, squared_euclidean_distance)) {
    return(rowSums((X1 - X2)^2))
  }
  if (identical(distance_function, euclidean_distance)) {
    return(sqrt(rowSums((X1 - X2)^2)))
  }
  if (identical(distance_function, manhattan_distance)) {
    return(rowSums(abs(X1 - X2)))
  }
  if (identical(distance_function, dot_product)) {
    return(rowSums(X1 * X2))
  }
  if (identical(distance_function, equality)) {
    return(as.numeric(rowSums(X1 == X2) == ncol(X1)))
  }
  # No known vectorised form: fall back to the generic, slower mapping.
  map_diag(distance_function, X1, X2)
}

#' Diagonal of the covariance matrix, vectorised when the kernel supports it
#'
#' @description
#' Default implementation: falls back to the always-correct, slower
#' [map_diag()] applied to [pairwise()]. Stationary and dot-product kernels
#' get a much faster method via [diag_distance()] + [shape()] -- the same
#' pattern as [pairwise_matrix()], one dimension lower.
#'
#' @param kernel A kernel object.
#' @param X1,X2 Numeric matrices with the same number of rows.
#' @return A numeric vector of length `nrow(X1)`.
#' @keywords internal
pairwise_diag <- function(kernel, X1, X2) {
  UseMethod("pairwise_diag")
}

#' @keywords internal
#' @exportS3Method
pairwise_diag.default <- function(kernel, X1, X2) {
  map_diag(function(x, y) pairwise(kernel, x, y), X1, X2)
}

#' @keywords internal
#' @exportS3Method
pairwise_diag.stationary_kernel <- function(kernel, X1, X2) {
  if (!has_s3_method("shape", kernel)) {
    return(pairwise_diag.default(kernel, X1, X2))
  }
  d <- diag_distance(kernel$distance_function, X1, X2)
  shape(kernel, d)
}

#' @keywords internal
#' @exportS3Method
pairwise_diag.dot_product_kernel <- function(kernel, X1, X2) {
  if (!has_s3_method("shape", kernel)) {
    return(pairwise_diag.default(kernel, X1, X2))
  }
  dp <- diag_distance(kernel$distance_function, X1, X2)
  shape(kernel, dp)
}
