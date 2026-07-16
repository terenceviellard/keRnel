#' Evaluate a kernel's covariance matrix
#'
#' @description
#' The top-level, user-facing entry point. Dispatches on the *class of the
#' kernel* (layer 4 of the architecture, see `dev/ameliorations/design/01-architecture.Rmd`):
#'
#' - The default method handles every base kernel: it coerces `x1`/`x2` to
#'   matrices, handles missing values, and delegates to [engine_eval()]
#'   (which in turn delegates to a mapping primitive and [pairwise()]).
#' - Composite/wrapper kernels (sums, products, exponentials, logarithms...)
#'   override `evaluate()` directly, bypassing `engine_eval()`/[pairwise()]
#'   entirely for themselves, and instead recombine the `evaluate()` results
#'   of their sub-kernels, forwarding `na_action`/`mask_variance` down via
#'   `...` so missing-value handling stays consistent however deeply a
#'   kernel is composed.
#'
#' Missing values (`NA`) in `x1`/`x2` are handled according to `na_action`:
#' - `"propagate"` (the default): `NA` in the input produces `NA` in the
#'   corresponding row/column of the result. The most predictable behaviour,
#'   but a covariance matrix containing `NA` cannot be inverted -- the
#'   caller must handle missing data explicitly before, e.g., a Cholesky
#'   factorisation.
#' - `"mask"`: a point with a missing value is treated as an independent
#'   observation of variance `mask_variance` (`1` by default), statistically
#'   neutral (adds no information about, and is not informed by, any other
#'   point) but keeps the resulting matrix positive-definite and therefore
#'   invertible. This is formalised here as an explicit `evaluate()` option
#'   rather than a separate engine, since engines should only encode *how* a
#'   result is vectorised, not domain-specific missing-data policy.
#'
#' With a kernel using [diagonal_engine()], the result is a vector (not a
#' matrix): `"propagate"` sets `NA` at index `i` whenever the aligned pair
#' `x1[i, ]`/`x2[i, ]` has a missing value; `"mask"` sets `mask_variance`
#' there instead.
#'
#' `evaluate()` also warns (not errors) if the matrix it is about to build
#' is large, since a dense `n1 x n2` matrix can consume significant memory
#' and compute time. Raise or lower `memory_warning_threshold`, or set it to
#' `Inf`, to adjust; this is a safety net, not a hard limit.
#'
#' A base kernel with a multi-valued hyperparameter (e.g.
#' `se_kernel(length_scale = c(1, 2, 3))`) errors if evaluated directly: such
#' a value is only meaningful sliced down by [batch_kernel()]/
#' [block_kernel()]/[input_specific_param_kernel()], which do so before
#' delegating to `evaluate()`. Evaluated directly, the vector would
#' otherwise be silently recycled across the covariance matrix by R's usual
#' recycling rules, producing a mathematically meaningless -- often
#' non-symmetric -- result with no warning at all.
#'
#' @param kernel A kernel object.
#' @param x1 A numeric vector or matrix of input points (see
#'   [as_matrix_input()]).
#' @param x2 An optional second set of input points; defaults to `x1`
#'   (the auto-covariance case).
#' @param na_action One of `"propagate"` (the default) or `"mask"`.
#' @param mask_variance A single non-negative numeric value, used when
#'   `na_action = "mask"`. Defaults to `1`.
#' @param memory_warning_threshold A single numeric value: `evaluate()` warns
#'   if `nrow(x1)` or `nrow(x2)` exceeds it. Defaults to `5000`; set to `Inf`
#'   to disable.
#' @param ... Passed on to methods (and, for composite kernels, forwarded to
#'   the `evaluate()` call of each sub-kernel). One internal argument
#'   travels this way: `.self_covariance`, the *self-covariance contract*.
#'   Whether the requested matrix is a self-covariance (`x2` missing or
#'   equal to `x1`) can only be decided at the level of the original call:
#'   wrappers that re-`evaluate()` their sub-kernel on transformed inputs
#'   (a single row of `x1`, one block's slice...) would mis-derive it
#'   locally, since a transformed `x1` no longer compares equal to `x2`.
#'   Such wrappers therefore decide it *once*, before transforming
#'   anything, and forward `.self_covariance` explicitly; the default
#'   method trusts a forwarded value over re-deriving its own. User code
#'   should never need to pass it.
#'
#' @return A numeric matrix of dimension `nrow(x1)` x `nrow(x2)` (or, with
#'   [diagonal_engine()], a numeric vector of length `nrow(x1)`).
#' @export
#' @examples
#' k <- se_kernel(length_scale = 1)
#' evaluate(k, c(0, 1, 2))
#' evaluate(k, c(0, 1), c(0, 10))
#'
#' # Keep the covariance matrix invertible despite missing data:
#' evaluate(k, c(0, NA, 2), na_action = "mask")
evaluate <- function(kernel, x1, x2 = NULL,
                      na_action = c("propagate", "mask"),
                      mask_variance = 1,
                      memory_warning_threshold = 5000,
                      ...) {
  na_action <- match.arg(na_action)

  n1 <- NROW(x1)
  n2 <- NROW(if (is.null(x2)) x1 else x2)
  if (n1 > memory_warning_threshold || n2 > memory_warning_threshold) {
    warning(
      "evaluate() is about to build a dense ", n1, "x", n2, " matrix (~",
      format(n1 * n2 * 8 / 1e6, digits = 3), " MB). ",
      "Consider diagonal_engine() if only the diagonal is needed, or raise ",
      "`memory_warning_threshold` (e.g. to Inf) to silence this warning.",
      call. = FALSE
    )
  }

  UseMethod("evaluate")
}

#' A hyperparameter with more than one value only has a well-defined meaning
#' when a batching/blocking wrapper mediates the evaluation -- `slice_kernel()`/
#' `slice_kernel_nd()`/`pair_slice_kernel()` (`R/wrappers.R`) always reduce it
#' to a single value (or a pair, for a [block_compatible][block_kernel()]
#' kernel) *before* the sliced kernel ever reaches `evaluate.default()`.
#' Reaching this function with a longer vector therefore means the kernel is
#' being evaluated directly, without such a wrapper: R would silently recycle
#' the vector across the covariance matrix (e.g. row-by-row, wrapping around),
#' producing a mathematically meaningless -- often non-symmetric -- result
#' with no error at all. Caught here, once, for every base kernel, rather
#' than duplicated in each `shape()`/`pairwise()` method.
#'
#' @keywords internal
.check_scalar_params <- function(kernel) {
  own <- Filter(function(field) inherits(field, "wrapped_param"), kernel)
  max_len <- if (isTRUE(attr(kernel, "block_compatible"))) 2L else 1L
  bad <- names(Filter(function(p) length(p$wrapped) > max_len, own))
  if (length(bad) > 0) {
    stop(
      "`", class(kernel)[1], "` is being evaluate()d directly, but hyperparameter(s) ",
      paste(bad, collapse = ", "), " have more than ", max_len, " value(s). A ",
      "multi-valued hyperparameter is only meaningful under batch_kernel(), ",
      "block_diag_kernel(), block_kernel(), or input_specific_param_kernel() ",
      "(which slice it down to a single value -- or a pair, for a ",
      "block_compatible kernel -- before evaluating); evaluated directly, it ",
      "would be silently recycled across the covariance matrix instead, ",
      "producing a mathematically meaningless (often non-symmetric) result. ",
      "Wrap this kernel in batch_kernel(), or give it a single shared value.",
      call. = FALSE
    )
  }
}

#' @export
evaluate.default <- function(kernel, x1, x2 = NULL,
                              na_action = "propagate",
                              mask_variance = 1,
                              memory_warning_threshold = 5000,
                              ...,
                              .self_covariance = NULL) {
  .check_scalar_params(kernel)

  # Trust a self-covariance status decided (and forwarded) by a wrapper
  # above us -- it saw the original, untransformed inputs; re-deriving it
  # here from already-transformed ones would get it wrong (see ?evaluate).
  is_self_covariance <- if (is.null(.self_covariance)) {
    is.null(x2) || identical(x1, x2)
  } else {
    isTRUE(.self_covariance)
  }
  if (is.null(x2)) {
    x2 <- x1
  }
  X1 <- as_matrix_input(x1)
  X2 <- as_matrix_input(x2)

  # `!is.finite()` catches NA/NaN/Inf/-Inf in one pass (is.finite() is FALSE
  # for all four) -- a row with an infinite coordinate is just as undefined
  # for a distance-based kernel as a missing one (e.g. `Inf - Inf` -> `NaN`
  # would otherwise leak into the covariance silently), so it gets the same
  # na_action treatment rather than a separate code path.
  na1 <- apply(X1, 1, function(row) any(!is.finite(row)))
  na2 <- apply(X2, 1, function(row) any(!is.finite(row)))
  X1[!is.finite(X1)] <- 0
  X2[!is.finite(X2)] <- 0

  result <- engine_eval(kernel$engine, kernel, X1, X2)

  if (na_action == "mask") {
    apply_na_mask(result, na1, na2, is_self_covariance, mask_variance)
  } else {
    apply_na_propagate(result, na1, na2)
  }
}

#' @keywords internal
apply_na_propagate <- function(result, na1, na2) {
  if (is.matrix(result)) {
    if (any(na1)) result[na1, ] <- NA
    if (any(na2)) result[, na2] <- NA
  } else {
    result[na1 | na2] <- NA
  }
  result
}

#' @keywords internal
apply_na_mask <- function(result, na1, na2, is_self_covariance, mask_variance) {
  if (is.matrix(result)) {
    if (any(na1)) result[na1, ] <- 0
    if (any(na2)) result[, na2] <- 0
    if (is_self_covariance) {
      diag_idx <- which(na1)
      result[cbind(diag_idx, diag_idx)] <- mask_variance
    }
  } else {
    result[na1 | na2] <- mask_variance
  }
  result
}

#' evaluate() with point identity propagated as dimnames
#'
#' @description
#' `evaluate()` itself never looks at `rownames(x1)`/`names(x1)` (a bare
#' matrix in, a bare matrix out) -- propagating point identity onto the
#' result is left to the caller. `evaluate_named()` is that propagation,
#' wrapped once: if `x1`/`x2` carry names (`rownames()` for a matrix or
#' data frame, `names()` for a bare vector), they become the result's
#' `rownames()`/`colnames()` (or, for a [batch_kernel()] array, the first
#' two dimnames; for a [diagonal_engine()] vector, its `names()`). With no
#' names on either side, the result is returned unchanged.
#'
#' @param kernel A kernel object.
#' @param x1 Inputs, as for [evaluate()]; a matrix/data frame's `rownames()`
#'   or a bare vector's `names()` trigger the propagation.
#' @param x2 Optional second set of inputs, as for [evaluate()].
#' @param ... Forwarded to [evaluate()].
#' @return The result of [evaluate()], with dimnames attached when available.
#' @export
#' @examples
#' x <- c(a = 0, b = 1, c = 2)
#' evaluate_named(se_kernel(length_scale = 1), x)
evaluate_named <- function(kernel, x1, x2 = NULL, ...) {
  result <- evaluate(kernel, x1, x2, ...)

  nm1 <- if (is.matrix(x1) || is.data.frame(x1)) rownames(x1) else names(x1)
  nm2 <- if (is.null(x2)) {
    nm1
  } else if (is.matrix(x2) || is.data.frame(x2)) {
    rownames(x2)
  } else {
    names(x2)
  }

  if (is.null(nm1) && is.null(nm2)) {
    return(result)
  }

  if (is.matrix(result)) {
    rownames(result) <- nm1
    colnames(result) <- nm2
  } else if (is.array(result)) {
    dn <- dimnames(result)
    if (is.null(dn)) dn <- vector("list", length(dim(result)))
    dn[[1]] <- nm1
    dn[[2]] <- nm2
    dimnames(result) <- dn
  } else {
    names(result) <- nm1
  }
  result
}
