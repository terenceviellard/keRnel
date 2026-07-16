#' Distance and similarity functions between input vectors
#'
#' @description
#' These functions implement the distance/comparison metrics used by
#' stationary and dot-product kernels (via their `distance_function`
#' argument). They all take two numeric vectors of the same length and
#' return a single scalar.
#'
#' @param x1,x2 Numeric vectors of the same length.
#' @param p Order of the norm for [minkowski_distance()] (`p >= 1`).
#' @param inv_cov_matrix Inverse covariance matrix for
#'   [mahalanobis_distance()].
#'
#' @return A single numeric value.
#' @name distances
NULL

#' @rdname distances
#' @export
#' @examples
#' euclidean_distance(c(0, 0), c(3, 4))
euclidean_distance <- function(x1, x2) sqrt(sum((x1 - x2)^2))

#' @rdname distances
#' @export
squared_euclidean_distance <- function(x1, x2) sum((x1 - x2)^2)

#' @rdname distances
#' @export
manhattan_distance <- function(x1, x2) sum(abs(x1 - x2))

#' @rdname distances
#' @export
chebyshev_distance <- function(x1, x2) max(abs(x1 - x2))

#' @rdname distances
#' @export
cosine_distance <- function(x1, x2) {
  1 - sum(x1 * x2) / (sqrt(sum(x1^2)) * sqrt(sum(x2^2)))
}

#' @rdname distances
#' @export
minkowski_distance <- function(x1, x2, p) {
  sum(abs(x1 - x2)^p)^(1 / p)
}

#' @rdname distances
#' @export
hamming_distance <- function(x1, x2) sum(x1 != x2)

#' @rdname distances
#' @export
jaccard_distance <- function(x1, x2) {
  intersection <- sum(pmin(x1, x2))
  union <- sum(pmax(x1, x2))
  1 - intersection / union
}

#' @rdname distances
#' @export
mahalanobis_distance <- function(x1, x2, inv_cov_matrix) {
  diff <- x1 - x2
  sqrt(as.numeric(t(diff) %*% inv_cov_matrix %*% diff))
}

#' @rdname distances
#' @export
canberra_distance <- function(x1, x2) {
  sum(abs(x1 - x2) / (abs(x1) + abs(x2) + 1e-10))
}

#' @rdname distances
#' @export
bray_curtis_distance <- function(x1, x2) {
  sum(abs(x1 - x2)) / sum(abs(x1 + x2) + 1e-10)
}

#' @rdname distances
#' @export
correlation_distance <- function(x1, x2) {
  x1_mean <- mean(x1)
  x2_mean <- mean(x2)
  numerator <- sum((x1 - x1_mean) * (x2 - x2_mean))
  denominator <- sqrt(sum((x1 - x1_mean)^2) * sum((x2 - x2_mean)^2))
  1 - numerator / (denominator + 1e-10)
}

#' @rdname distances
#' @export
dot_product <- function(x1, x2) sum(x1 * x2)

#' @rdname distances
#' @export
#' @examples
#' equality(c(1, 2), c(1, 2))
#' equality(c(1, 2), c(1, 3))
equality <- function(x1, x2) as.numeric(isTRUE(all(x1 == x2)))

# Registry used by resolve_distance_function(). Only the simple two-argument
# distances are listed here -- minkowski_distance()/mahalanobis_distance()
# need an extra argument (`p`/`inv_cov_matrix`), so a bare string shortcut
# would not be enough to fully specify them.
.distance_registry <- list(
  euclidean = euclidean_distance,
  squared_euclidean = squared_euclidean_distance,
  manhattan = manhattan_distance,
  chebyshev = chebyshev_distance,
  cosine = cosine_distance,
  hamming = hamming_distance,
  jaccard = jaccard_distance,
  canberra = canberra_distance,
  bray_curtis = bray_curtis_distance,
  correlation = correlation_distance,
  dot_product = dot_product,
  equality = equality
)

#' Resolve a distance function from a function or a shortcut name
#'
#' @description
#' Lets kernel constructors accept either a distance function directly, or a
#' convenience string shortcut (e.g. `"euclidean"`), resolved against an
#' internal registry. Exported so a custom `stationary_kernel`/
#' `dot_product_kernel` constructor (see `vignette("b-intermediate")`,
#' "Writing your own kernel") can support the same shortcut for its own
#' `distance_function` argument.
#'
#' This shortcut is a deliberate ergonomic addition for R, mirroring
#' idioms such as `stats::dist(method = "euclidean")`.
#'
#' @param x A function, or a single string naming one of the distances in
#'   [distances] (run `names(keRnel:::.distance_registry)` for the full
#'   list).
#' @return A distance function taking `(x1, x2)`.
#' @export
#' @examples
#' resolve_distance_function("manhattan")
resolve_distance_function <- function(x) {
  if (is.function(x)) {
    return(x)
  }
  if (is.character(x) && length(x) == 1 && x %in% names(.distance_registry)) {
    return(.distance_registry[[x]])
  }
  stop(
    "`distance_function` must be a function or one of: ",
    paste(names(.distance_registry), collapse = ", "),
    ".",
    call. = FALSE
  )
}
