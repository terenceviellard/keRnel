#' Pairwise kernel formula
#'
#' @description
#' Computes the covariance between a single pair of points. This is the
#' *only* thing a new base kernel needs to implement (layer 1 of the
#' architecture, see `dev/ameliorations/design/01-architecture.Rmd`) -- vectorisation over
#' many points is handled separately by the mapping/engine layers, via
#' [evaluate()].
#'
#' @param kernel A kernel object.
#' @param x1,x2 Numeric vectors representing a single point each (one row of
#'   the matrices passed to [evaluate()]).
#'
#' @return A single numeric value.
#' @export
pairwise <- function(kernel, x1, x2) {
  UseMethod("pairwise")
}
