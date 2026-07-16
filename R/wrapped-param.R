#' Wrapped hyperparameter
#'
#' @description
#' The atomic unit of an optimisable kernel hyperparameter. It bundles
#' together:
#' - `wrapped`: the value in the unconstrained ("wrapped") space.
#' - `parametrisation`: the [parametrisation] object used to `wrap()`/
#'   `unwrap()` it.
#' - `frozen`: whether this parameter should be excluded from optimisation
#'   (see [freeze()] and [get_trainable_params()]).
#' - `sliceable`: whether a multi-valued vector for this parameter indexes a
#'   batch/block axis (and may therefore be sliced by [batch_kernel()] and
#'   friends, see [slice_kernel()]), or something else entirely -- e.g.
#'   [ard_kernel()]'s `length_scales`, whose entries index *input
#'   dimensions*, never batch elements. Non-sliceable parameters are ignored
#'   by batch-axis auto-detection and rejected if named explicitly in a
#'   `batch_axes`-style argument.
#'
#' `wrapped_param` bundles a hyperparameter's value and its parametrisation
#' into a single, S3-dispatchable object, rather than keeping them as
#' separate fields linked only by a naming convention. This makes generic
#' tree-walking functions (`get_params()`, `kupdate()`, `slice_kernel()`...)
#' robust: they dispatch on `inherits(x, "wrapped_param")` rather than
#' guessing field names.
#'
#' The *intent* behind `frozen` -- a parameter that participates in the
#' computation normally but must never be optimised -- is fully captured by
#' `frozen = TRUE` here, with no need for a dedicated parametrisation class.
#' `sliceable` follows the same philosophy for the batching axis: metadata
#' carried directly by the parameter itself, no parallel structure to keep
#' in sync.
#'
#' @param value A numeric value or vector, in the constrained (natural) space.
#' @param parametrisation A [parametrisation] object. Defaults to
#'   [log_exp_parametrisation()] (the most common case: a strictly positive
#'   hyperparameter).
#' @param frozen Logical, whether the parameter is excluded from optimisation.
#' @param sliceable Logical, whether a multi-valued vector for this parameter
#'   may be interpreted as one-value-per-batch-element. Defaults to `TRUE`.
#'
#' @return An object of class `wrapped_param`.
#' @export
#' @examples
#' p <- wrapped_param(2.5)
#' unwrap_param(p)
#'
#' p_frozen <- wrapped_param(1, frozen = TRUE)
#' p_frozen$frozen
wrapped_param <- function(value,
                           parametrisation = log_exp_parametrisation(),
                           frozen = FALSE,
                           sliceable = TRUE) {
  structure(
    list(
      wrapped = wrap(parametrisation, value),
      parametrisation = parametrisation,
      frozen = isTRUE(frozen),
      sliceable = !isFALSE(sliceable)
    ),
    class = "wrapped_param"
  )
}

#' Read the constrained value of a wrapped parameter
#'
#' @param x An object of class `wrapped_param`.
#' @return The value in the natural (constrained) space.
#' @export
#' @examples
#' unwrap_param(wrapped_param(2.5))
unwrap_param <- function(x) {
  stopifnot(inherits(x, "wrapped_param"))
  unwrap(x$parametrisation, x$wrapped)
}

#' @export
print.wrapped_param <- function(x, ...) {
  cat(format(x), "\n", sep = "")
  invisible(x)
}

#' @export
format.wrapped_param <- function(x, ...) {
  value <- unwrap_param(x)
  suffix <- if (isTRUE(x$frozen)) " [frozen]" else ""
  # Always collapse to a single string, even for a vector-valued parameter
  # (e.g. ard_kernel()'s length_scales, or a batch_kernel()-sliceable
  # per-element value) -- format(<vector>) alone would return one string
  # per element, breaking callers (format.kernel()) that expect exactly one
  # string per hyperparameter.
  body <- if (length(value) > 1) {
    paste0("[", paste(format(value, digits = 4), collapse = ", "), "]")
  } else {
    format(value, digits = 4)
  }
  paste0(body, suffix)
}
