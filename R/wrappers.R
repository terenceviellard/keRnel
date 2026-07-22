# Wrapper kernels: unary compositions, stored as $subkernels (length 1),
# with no engine of their own -- they override evaluate() directly, applying
# a transform to their sub-kernel's evaluate() result (layer 4 of the
# architecture, see dev/ameliorations/design/01-architecture.Rmd and R/operators.R).
#
# active_dims_kernel() and ard_kernel() transform the *inputs* before
# delegating to the inner kernel's evaluate(). batch_kernel(),
# input_specific_param_kernel(), block_diag_kernel(), and block_kernel()
# instead evaluate() the inner kernel once per batch/block/input element
# (via slice_kernel()/pair_slice_kernel()) and reassemble the results.

#' Exponential of a kernel
#'
#' @description
#' \deqn{k(x_1, x_2) = \exp(k_{\mathrm{inner}}(x_1, x_2))}
#'
#' @param inner A kernel object, or a numeric value (coerced via
#'   [as_kernel()]).
#'
#' @return An object of class `exp_kernel`/`wrapper_kernel`/`kernel`.
#' @export
#' @examples
#' k <- exp_kernel(se_kernel(length_scale = 1))
#' evaluate(k, c(0, 1, 2))
exp_kernel <- function(inner) {
  kernel <- structure(
    list(subkernels = list(as_kernel(inner))),
    class = c("exp_kernel", "wrapper_kernel", "kernel")
  )
  validate_kernel(kernel)
  kernel
}

#' @export
evaluate.exp_kernel <- function(kernel, x1, x2 = NULL, ...) {
  exp(evaluate(kernel$subkernels[[1]], x1, x2, ...))
}

# d(exp(k_inner))/dtheta = exp(k_inner) * dk_inner/dtheta -- the ordinary
# chain rule, reusing .grad_times() (R/kernel-grad.R) so a vector-valued
# hyperparameter's list-of-matrices entry (e.g. feature_kernel()'s
# length_scale) is handled the same way as a plain matrix one.

#' @keywords internal
#' @exportS3Method
kernel_grad_raw.exp_kernel <- function(kernel, x1, x2) {
  inner <- kernel$subkernels[[1]]
  k_outer <- exp(evaluate(inner, x1, x2))
  lapply(kernel_grad_raw(inner, x1, x2), .grad_times, b = k_outer)
}

#' @export
format.exp_kernel <- function(x, ...) {
  paste0("Exp(", format(x$subkernels[[1]]), ")")
}

#' Logarithm of a kernel
#'
#' @description
#' \deqn{k(x_1, x_2) = \log(k_{\mathrm{inner}}(x_1, x_2))}
#'
#' Note: `k_inner`'s output must be strictly positive everywhere for the
#' result to be defined; this is not checked at construction (it generally
#' cannot be, since it depends on the inputs `evaluate()` is later called
#' with), only at evaluation time (`NaN`, with an R warning, as usual for
#' `log()` of a non-positive value).
#'
#' @param inner A kernel object, or a numeric value (coerced via
#'   [as_kernel()]).
#'
#' @return An object of class `log_kernel`/`wrapper_kernel`/`kernel`.
#' @export
#' @examples
#' k <- log_kernel(variance_kernel(variance = 2))
#' evaluate(k, c(0, 1, 2))
log_kernel <- function(inner) {
  kernel <- structure(
    list(subkernels = list(as_kernel(inner))),
    class = c("log_kernel", "wrapper_kernel", "kernel")
  )
  validate_kernel(kernel)
  kernel
}

#' @export
evaluate.log_kernel <- function(kernel, x1, x2 = NULL, ...) {
  log(evaluate(kernel$subkernels[[1]], x1, x2, ...))
}

# d(log(k_inner))/dtheta = dk_inner/dtheta / k_inner.

#' @keywords internal
#' @exportS3Method
kernel_grad_raw.log_kernel <- function(kernel, x1, x2) {
  inner <- kernel$subkernels[[1]]
  k_inner <- evaluate(inner, x1, x2)
  lapply(kernel_grad_raw(inner, x1, x2), function(g) .grad_times(g, 1 / k_inner))
}

#' @export
format.log_kernel <- function(x, ...) {
  paste0("Log(", format(x$subkernels[[1]]), ")")
}

#' Select input dimensions before applying a kernel
#'
#' @description
#' Restricts the inner kernel to a subset of input dimensions (columns),
#' selected before `evaluate()` is delegated to it:
#'
#' \deqn{k(x_1, x_2) = k_{\mathrm{inner}}(x_1[, \mathrm{active\_dims}], x_2[, \mathrm{active\_dims}])}
#'
#' Note: this wrapper *must* be the outer-most one (it should not itself be
#' nested inside another wrapper) -- a documented usage contract, not
#' automatically verified (verifying it in general would require inspecting
#' whether a kernel is nested inside anything else, which is not knowable
#' from the kernel object itself).
#'
#' @param inner A kernel object, or a numeric value (coerced via
#'   [as_kernel()]).
#' @param active_dims A vector of positive integers: the (1-based) column
#'   indices of the input dimensions to keep.
#'
#' @return An object of class `active_dims_kernel`/`wrapper_kernel`/`kernel`.
#' @export
#' @examples
#' k <- active_dims_kernel(se_kernel(length_scale = 1), active_dims = 2)
#' X <- matrix(c(0, 10, 1, 11, 2, 12), ncol = 2, byrow = TRUE)
#' evaluate(k, X)  # only column 2 (values 10, 11, 12) is seen by se_kernel
active_dims_kernel <- function(inner, active_dims) {
  active_dims <- as.integer(active_dims)
  if (length(active_dims) == 0 || any(active_dims < 1)) {
    stop("`active_dims` must be one or more positive integers (1-based column indices).", call. = FALSE)
  }

  kernel <- structure(
    list(subkernels = list(as_kernel(inner)), active_dims = active_dims),
    class = c("active_dims_kernel", "wrapper_kernel", "kernel")
  )
  validate_kernel(kernel)
  kernel
}

#' @export
evaluate.active_dims_kernel <- function(kernel, x1, x2 = NULL, ...) {
  X1 <- as_matrix_input(x1)
  if (max(kernel$active_dims) > ncol(X1)) {
    stop(
      "`active_dims` (", paste(kernel$active_dims, collapse = ", "), ") is out of bounds ",
      "for an input with ", ncol(X1), " column(s).",
      call. = FALSE
    )
  }
  X1 <- X1[, kernel$active_dims, drop = FALSE]
  X2 <- if (is.null(x2)) NULL else as_matrix_input(x2)[, kernel$active_dims, drop = FALSE]
  evaluate(kernel$subkernels[[1]], X1, X2, ...)
}

# active_dims_kernel() has no hyperparameter of its own (active_dims is a
# fixed integer index, not optimisable) -- it just selects columns before
# delegating, so its gradient is the inner kernel's gradient on the same
# selected columns, no new derivative math needed.

#' @keywords internal
#' @exportS3Method
kernel_grad_raw.active_dims_kernel <- function(kernel, x1, x2) {
  X1 <- as_matrix_input(x1)
  X2 <- as_matrix_input(if (is.null(x2)) x1 else x2)
  kernel_grad_raw(
    kernel$subkernels[[1]],
    X1[, kernel$active_dims, drop = FALSE],
    X2[, kernel$active_dims, drop = FALSE]
  )
}

#' @export
format.active_dims_kernel <- function(x, ...) {
  paste0(format(x$subkernels[[1]]), "[dims=", paste(x$active_dims, collapse = ","), "]")
}

#' Automatic Relevance Determination (ARD)
#'
#' @description
#' Wraps a stationary kernel constructor to give every input dimension its
#' own `length_scale`:
#'
#' \deqn{k(x_1, x_2) = k_{\mathrm{inner}}(x_1 / \mathrm{length\_scales}, x_2 / \mathrm{length\_scales})}
#'
#' `ard_kernel()` takes a kernel **constructor** (a function), not a
#' pre-built instance. It builds the inner kernel itself with
#' `length_scale = 1` and immediately [freeze()]s it -- the contract that
#' every per-dimension `length_scale` must start at `1` and stay frozen
#' cannot be violated by the caller, because the caller never gets to build
#' (or therefore mis-build) the inner kernel.
#'
#' (Residual risk, accepted and documented: [kupdate()] does not itself
#' check `frozen` before overwriting a value, so
#' `kupdate(ard, length_scale = 5)` would still bypass this protection if
#' called directly on the inner kernel's hyperparameter name -- a
#' deliberate, explicitly documented limitation rather than an oversight.)
#'
#' @param constructor A kernel constructor function, e.g. [se_kernel()],
#'   taking a `length_scale` argument (or whichever argument name `active`
#'   specifies).
#' @param length_scales A vector of strictly positive numeric values, one
#'   per input dimension.
#' @param ... Additional arguments passed to `constructor` (e.g. `period`
#'   for [periodic_kernel()]).
#' @param active The name of `constructor`'s length-scale argument. Defaults
#'   to `"length_scale"`.
#' @param length_scales_parametrisation A [parametrisation] for
#'   `length_scales`. Defaults to [log_exp_parametrisation()].
#'
#' @return An object of class `ard_kernel`/`wrapper_kernel`/`kernel`.
#' @export
#' @examples
#' k <- ard_kernel(se_kernel, length_scales = c(1, 2, 3))
#' X <- matrix(rnorm(9), ncol = 3)
#' evaluate(k, X)
#'
#' # Additional constructor arguments are passed through `...`:
#' k2 <- ard_kernel(periodic_kernel, length_scales = c(1, 2), period = 5)
ard_kernel <- function(constructor,
                        length_scales,
                        ...,
                        active = "length_scale",
                        length_scales_parametrisation = log_exp_parametrisation()) {
  constructor_args <- list(...)
  constructor_args[[active]] <- 1
  inner <- do.call(constructor, constructor_args)
  inner <- freeze(inner, active)

  kernel <- structure(
    list(
      subkernels = list(inner),
      # sliceable = FALSE: length_scales has one value per *input
      # dimension*, not per batch element -- batch_kernel() and friends
      # must never mistake it for a batch axis, even when its length
      # happens to coincide with a batch_size.
      length_scales = wrapped_param(length_scales, length_scales_parametrisation,
                                     sliceable = FALSE)
    ),
    class = c("ard_kernel", "wrapper_kernel", "kernel")
  )
  validate_kernel(kernel)
  kernel
}

#' @export
evaluate.ard_kernel <- function(kernel, x1, x2 = NULL, ...) {
  length_scales <- unwrap_param(kernel$length_scales)
  X1 <- as_matrix_input(x1)
  if (ncol(X1) != length(length_scales)) {
    stop(
      "Input has ", ncol(X1), " column(s) but `length_scales` has length ",
      length(length_scales), ".",
      call. = FALSE
    )
  }
  X1 <- sweep(X1, 2, length_scales, "/")
  X2 <- if (is.null(x2)) NULL else sweep(as_matrix_input(x2), 2, length_scales, "/")
  evaluate(kernel$subkernels[[1]], X1, X2, ...)
}

# ard_kernel() feeds its inner kernel distance_function(x1 / length_scales,
# x2 / length_scales) -- length_scales therefore enters through the
# *distance*, not directly as one of the inner kernel's own hyperparameters
# (which shape_grad() differentiates against). d(distance)/d(length_scales[k])
# is derived here via the chain rule through u = x / length_scales,
# dimension by dimension, for the vectorised distance functions
# distance_matrix() itself knows a closed form for (besides equality(), a
# non-differentiable step function): writing f(u1, u2) for the distance in
# rescaled space,
#   d(f)/d(length_scales[k]) = df/du1[k] * du1[k]/d(length_scales[k])
#                             + df/du2[k] * du2[k]/d(length_scales[k])
# with du1[k]/d(length_scales[k]) = -x1[k] / length_scales[k]^2 (and
# likewise for u2) -- only dimension k's own length_scale affects u1[k]/
# u2[k], so no cross-dimension terms appear.

#' @keywords internal
.ard_distance_grad <- function(distance_function, X1, X2, length_scales) {
  is_squared_euclidean <- identical(distance_function, squared_euclidean_distance)
  is_euclidean <- identical(distance_function, euclidean_distance)
  is_manhattan <- identical(distance_function, manhattan_distance)
  is_dot_product <- identical(distance_function, dot_product)
  if (!any(is_squared_euclidean, is_euclidean, is_manhattan, is_dot_product)) {
    stop(
      "ard_kernel() has no analytic gradient for this `distance_function` -- only ",
      "squared_euclidean_distance(), euclidean_distance(), manhattan_distance(), ",
      "and dot_product() are supported.",
      call. = FALSE
    )
  }

  U1 <- sweep(X1, 2, length_scales, "/")
  U2 <- sweep(X2, 2, length_scales, "/")
  d <- if (is_euclidean) distance_matrix(distance_function, U1, U2) else NULL

  lapply(seq_len(ncol(X1)), function(k) {
    diff_k <- outer(U1[, k], U2[, k], "-")
    if (is_squared_euclidean) {
      df_du1 <- 2 * diff_k
      df_du2 <- -df_du1
    } else if (is_euclidean) {
      df_du1 <- ifelse(d == 0, 0, diff_k / d)
      df_du2 <- -df_du1
    } else if (is_manhattan) {
      df_du1 <- sign(diff_k)
      df_du2 <- -df_du1
    } else {
      df_du1 <- matrix(U2[, k], nrow = nrow(U1), ncol = nrow(U2), byrow = TRUE)
      df_du2 <- matrix(U1[, k], nrow = nrow(U1), ncol = nrow(U2))
    }
    du1_k <- -X1[, k] / length_scales[k]^2
    du2_k <- -X2[, k] / length_scales[k]^2
    sweep(df_du1, 1, du1_k, "*") + sweep(df_du2, 2, du2_k, "*")
  })
}

#' @keywords internal
#' @exportS3Method
kernel_grad_raw.ard_kernel <- function(kernel, x1, x2) {
  inner <- kernel$subkernels[[1]]
  length_scales <- unwrap_param(kernel$length_scales)
  X1 <- as_matrix_input(x1)
  X2 <- as_matrix_input(if (is.null(x2)) x1 else x2)
  X1s <- sweep(X1, 2, length_scales, "/")
  X2s <- sweep(X2, 2, length_scales, "/")

  own_grad <- list()
  if (!isTRUE(kernel$length_scales$frozen)) {
    if (!has_s3_method("shape", inner) || !has_s3_method("shape_grad_d", inner)) {
      stop(
        "kernel_grad() is not implemented for ard_kernel(", class(inner)[1], ", ...) -- ",
        "analytic ARD gradients require the wrapped kernel to have a shape()/",
        "shape_grad_d() method (the base stationary/dot-product kernels do; ",
        "feature_kernel() and composite kernels do not).",
        call. = FALSE
      )
    }
    d <- distance_matrix(inner$distance_function, X1s, X2s)
    dd_dls <- .ard_distance_grad(inner$distance_function, X1, X2, length_scales)
    dshape_dd <- shape_grad_d(inner, d)
    own_grad <- list(length_scales = lapply(dd_dls, function(g) dshape_dd * g))
  }

  c(own_grad, kernel_grad_raw(inner, X1s, X2s))
}

#' @export
format.ard_kernel <- function(x, ...) {
  paste0("ARD(", format(x$subkernels[[1]]), ", length_scales = ", format(x$length_scales), ")")
}

# batch_kernel()/input_specific_param_kernel() both need to pick out, for a
# single batch/input index, the "slice" of a kernel whose hyperparameters
# vary per index. Which hyperparameters vary along the wrapper's axis is an
# explicit *contract*, resolved once at construction by .resolve_axes()
# below (either declared by name via `batch_axes`-style arguments, or
# auto-detected from vector lengths -- with strict validation either way),
# then stored on the wrapper and re-checked at evaluate() time (so a later
# kupdate() cannot silently invalidate it). Slicing itself
# (slice_kernel()/pair_slice_kernel()) only ever touches the resolved axis
# names -- never "any multi-valued parameter", which would confuse a
# per-dimension vector (e.g. ard_kernel()'s length_scales) with a per-batch
# one. Frozen/sliceable flags live directly on each wrapped_param rather
# than in a separate parallel structure, so nested batching composes
# without special-casing.

#' Resolve which hyperparameters vary along a wrapper's batch/block axis
#'
#' @description
#' The single place where the "which hyperparameter is batched?" question is
#' answered, shared by [batch_kernel()], [block_diag_kernel()],
#' [block_kernel()] and [input_specific_param_kernel()]. Returns the axis
#' names as a character vector, after validating them against the kernel
#' tree; called both at construction time and again at `evaluate()` time
#' (a [kupdate()] in between may have changed a value's length).
#'
#' Two modes:
#' - `axes = NULL` (auto-detection): every *sliceable* multi-valued
#'   `wrapped_param` is taken to vary along this wrapper's axis, and must
#'   therefore have exactly `n` values -- any other length is an error
#'   (silently recycling or truncating a hyperparameter vector is how
#'   covariance matrices go quietly wrong). Parameters marked
#'   `sliceable = FALSE` (e.g. [ard_kernel()]'s per-*dimension*
#'   `length_scales`) are skipped: their length indexes something else
#'   entirely.
#' - `axes = <character vector>`: only the named hyperparameters vary along
#'   this axis; every occurrence must have length `n` (or 1, i.e. shared).
#'   Other multi-valued parameters are left alone -- they belong to some
#'   *other* wrapper's axis (this is what makes nested batching with
#'   different sizes per level work).
#'
#' @param inner A kernel object (the wrapper's sub-kernel).
#' @param n A positive integer: the wrapper's axis size.
#' @param axes `NULL`, or a character vector of hyperparameter names.
#' @param axes_arg,size_arg The wrapper's argument names, for error messages
#'   (e.g. `"batch_axes"`, `"batch_size"`).
#' @return A character vector of resolved axis names (possibly empty).
#' @keywords internal
.resolve_axes <- function(inner, n, axes, axes_arg, size_arg) {
  params <- flatten_params(inner)
  lens <- vapply(params, function(p) length(p$wrapped), integer(1))
  sliceable <- vapply(params, function(p) !isFALSE(p$sliceable), logical(1))

  if (is.null(axes)) {
    multi <- lens > 1 & sliceable
    bad <- unique(names(params)[multi & lens != n])
    if (length(bad) > 0) {
      bad_lens <- lens[names(params) %in% bad & multi]
      stop(
        "Hyperparameter(s) ", paste(bad, collapse = ", "), " have ",
        paste(unique(bad_lens), collapse = "/"), " value(s) but `", size_arg,
        "` is ", n, ". Give each per-element hyperparameter exactly ", n,
        " values, or declare which hyperparameters vary along this axis ",
        "explicitly via `", axes_arg, "` (multi-valued hyperparameters not ",
        "named there are left to other wrappers, e.g. a nested batch_kernel()).",
        call. = FALSE
      )
    }
    return(unique(names(params)[multi & lens == n]))
  }

  axes <- unique(as.character(axes))
  if (length(axes) == 0) {
    return(character(0))
  }

  available <- unique(names(params))
  missing <- setdiff(axes, available)
  if (length(missing) > 0) {
    stop(
      "No such hyperparameter(s) in this kernel: ", paste(missing, collapse = ", "),
      ". Available: ", paste(available, collapse = ", "), ".",
      call. = FALSE
    )
  }

  not_sliceable <- intersect(axes, unique(names(params)[!sliceable]))
  if (length(not_sliceable) > 0) {
    stop(
      "Hyperparameter(s) ", paste(not_sliceable, collapse = ", "), " cannot be ",
      "used as a batch/block axis: their values index something other than ",
      "batch elements (e.g. ard_kernel()'s `length_scales` has one value per ",
      "*input dimension*).",
      call. = FALSE
    )
  }

  for (a in axes) {
    occ_lens <- lens[names(params) == a]
    bad <- unique(occ_lens[!(occ_lens %in% c(1L, n))])
    if (length(bad) > 0) {
      stop(
        "Hyperparameter `", a, "` is declared in `", axes_arg, "` but has ",
        paste(bad, collapse = "/"), " value(s) where `", size_arg, "` is ", n,
        ". A hyperparameter varying along this axis needs exactly ", n,
        " values (or a single shared one).",
        call. = FALSE
      )
    }
  }
  axes
}

# Multi-axis batching: a single wrapper handling several *named* batch
# dimensions at once (e.g. `c(K = 3, O = 2)`), so that one hyperparameter can
# vary along the full Cartesian product of two or more axes simultaneously --
# something a chain of nested batch_kernel() calls cannot express, because
# each level's slice_kernel() collapses an axis hyperparameter to a scalar
# before the inner level ever sees it (see dev/ameliorations/design/07-contrat-multi-axes.Rmd).
# `batch_size` becomes a named integer vector, and `batch_axes` becomes a
# named *list* (hyperparameter name -> character vector of the axis names,
# a subset of `names(batch_size)`, that hyperparameter varies along) instead
# of a plain character vector. Auto-detection is deliberately not supported
# here: with several axes in play, a hyperparameter's length alone is not
# enough to tell which axis (or which combination) it belongs to.
#
# Flattening convention (documented in batch_kernel()'s roxygen): a
# hyperparameter declared to vary along axes `c(a1, a2, ...)` stores its
# values in column-major order over `batch_size[c(a1, a2, ...)]` -- exactly
# the order `array(seq_len(prod(sizes)), dim = sizes)` would fill, i.e. `a1`
# is the fastest-varying axis. The wrapper's own output array follows the
# same convention over the *full* `names(batch_size)` order.

#' @keywords internal
.validate_multi_batch_size <- function(batch_size, size_arg) {
  nms <- names(batch_size)
  if (is.null(nms) || any(nms == "") || anyDuplicated(nms)) {
    stop(
      "A multi-axis `", size_arg, "` (length > 1) must be a named integer ",
      "vector with unique, non-empty names, e.g. `c(K = 3, O = 2)`.",
      call. = FALSE
    )
  }
  sizes <- suppressWarnings(as.integer(batch_size))
  if (anyNA(sizes) || any(sizes < 1)) {
    stop("Every axis size in `", size_arg, "` must be a positive integer.", call. = FALSE)
  }
  names(sizes) <- nms
  sizes
}

#' Resolve which hyperparameters vary along which axis of a multi-axis batch wrapper
#'
#' @description
#' The multi-axis counterpart of [.resolve_axes()]. `batch_axes` must be an
#' explicit named list (hyperparameter name -> character vector of axis
#' names, a subset of `names(batch_size)`) -- auto-detection is not
#' supported, since a hyperparameter's length alone cannot disambiguate
#' which axis, or which combination of axes, it belongs to once more than
#' one axis is in play. Re-validated at every `evaluate()` call, exactly
#' like [.resolve_axes()].
#'
#' @param inner A kernel object (the wrapper's sub-kernel).
#' @param batch_size A named integer vector (axis name -> axis size).
#' @param batch_axes A named list (hyperparameter name -> character vector
#'   of axis names).
#' @param axes_arg,size_arg The wrapper's argument names, for error messages.
#' @return `batch_axes`, unchanged, if every declaration is valid.
#' @keywords internal
.resolve_axes_multi <- function(inner, batch_size, batch_axes, axes_arg, size_arg) {
  if (is.null(batch_axes)) {
    stop(
      "Multi-axis `", size_arg, "` (", paste(names(batch_size), collapse = ", "),
      ") requires an explicit `", axes_arg, "`: a named list mapping each ",
      "varying hyperparameter to the subset of axis names it varies along. ",
      "Auto-detection is not supported for multi-axis batching -- see ?batch_kernel.",
      call. = FALSE
    )
  }
  if (!is.list(batch_axes) ||
        (length(batch_axes) > 0 &&
           (is.null(names(batch_axes)) || any(names(batch_axes) == "")))) {
    stop(
      "`", axes_arg, "` must be a named list (hyperparameter name -> ",
      "character vector of axis names) when `", size_arg, "` has more than ",
      "one axis. An empty list (`list()`) is valid: no hyperparameter ",
      "varies, every one is shared across the whole grid of axes.",
      call. = FALSE
    )
  }
  dup_keys <- unique(names(batch_axes)[duplicated(names(batch_axes))])
  if (length(dup_keys) > 0) {
    stop(
      "`", axes_arg, "` has more than one entry for hyperparameter(s): ",
      paste(dup_keys, collapse = ", "), ". `[[` lookup by name would silently ",
      "keep only the first and ignore the rest -- merge them into a single ",
      "entry instead.",
      call. = FALSE
    )
  }

  params <- flatten_params(inner)
  available <- unique(names(params))
  resolved <- batch_axes

  for (pname in names(batch_axes)) {
    axes_h <- unique(as.character(batch_axes[[pname]]))
    if (length(axes_h) == 0) {
      stop(
        "`", axes_arg, "[[\"", pname, "\"]]` is empty. To leave `", pname,
        "` fully shared across every axis, omit it from `", axes_arg,
        "` entirely rather than naming it with zero axes.",
        call. = FALSE
      )
    }

    unknown <- setdiff(axes_h, names(batch_size))
    if (length(unknown) > 0) {
      stop(
        "`", axes_arg, "[[\"", pname, "\"]]` names unknown axis/axes: ",
        paste(unknown, collapse = ", "), ". Available axes: ",
        paste(names(batch_size), collapse = ", "), ".",
        call. = FALSE
      )
    }
    if (!(pname %in% available)) {
      stop(
        "No such hyperparameter in this kernel: ", pname, ". Available: ",
        paste(available, collapse = ", "), ".",
        call. = FALSE
      )
    }

    occ <- params[names(params) == pname]
    sliceable <- vapply(occ, function(p) !isFALSE(p$sliceable), logical(1))
    if (any(!sliceable)) {
      stop(
        "Hyperparameter `", pname, "` cannot be used as a batch axis: its ",
        "values index something other than batch elements (e.g. ",
        "ard_kernel()'s `length_scales`).",
        call. = FALSE
      )
    }

    expected <- prod(batch_size[axes_h])
    lens <- vapply(occ, function(p) length(p$wrapped), integer(1))
    bad <- unique(lens[!(lens %in% c(1L, expected))])
    if (length(bad) > 0) {
      stop(
        "Hyperparameter `", pname, "` is declared to vary along axes ",
        paste(axes_h, collapse = " x "), " (", expected, " combination(s)) ",
        "but has ", paste(bad, collapse = "/"), " value(s). Give it exactly ",
        expected, " values (column-major order over ",
        paste0("batch_size[c(", paste(sprintf('"%s"', axes_h), collapse = ", "), ")]"),
        "), or a single shared one.",
        call. = FALSE
      )
    }

    # Store the deduplicated axis vector (not the raw, possibly-repeated
    # input) -- slice_kernel_nd() indexes by this contract directly, and a
    # repeated axis name there would double-weight that axis in the
    # column-major flat-index arithmetic, reading past the hyperparameter's
    # actual length.
    resolved[[pname]] <- axes_h
  }

  resolved
}

#' Extract one Cartesian combination of a multi-axis batch wrapper's hyperparameters
#'
#' @description
#' The multi-axis counterpart of [slice_kernel()]. `combo` gives one index
#' per wrapper axis (a named integer vector over `names(batch_size)`); each
#' hyperparameter named in `axis_spec` is replaced by the element of its
#' `$wrapped` vector at the flattened (column-major) index implied by
#' `combo` restricted to that hyperparameter's own declared axes -- see the
#' flattening convention documented in [batch_kernel()].
#'
#' @param kernel A kernel object.
#' @param combo A named integer vector (axis name -> current 1-based index).
#' @param axis_spec A named list (hyperparameter name -> character vector of
#'   axis names), the contract resolved by `.resolve_axes_multi()`.
#' @param batch_size A named integer vector (axis name -> axis size).
#' @return A new kernel object, with every axis hyperparameter sliced at the
#'   combination `combo` implies.
#' @keywords internal
slice_kernel_nd <- function(kernel, combo, axis_spec, batch_size) {
  for (name in names(kernel)) {
    field <- kernel[[name]]
    if (inherits(field, "wrapped_param") && name %in% names(axis_spec) &&
          length(field$wrapped) > 1) {
      axes_h <- axis_spec[[name]]
      sizes_h <- as.integer(batch_size[axes_h])
      idx_h <- as.integer(combo[axes_h])
      flat_idx <- if (length(sizes_h) == 1) {
        idx_h
      } else {
        1L + sum((idx_h - 1L) * c(1L, cumprod(sizes_h))[seq_along(sizes_h)])
      }
      kernel[[name]]$wrapped <- field$wrapped[flat_idx]
    }
  }
  if (!is.null(kernel$subkernels)) {
    kernel$subkernels <- lapply(kernel$subkernels, slice_kernel_nd,
                                 combo = combo, axis_spec = axis_spec,
                                 batch_size = batch_size)
  }
  kernel
}

#' Extract one batch/input element of a kernel's hyperparameters
#'
#' @description
#' Walks `kernel` (recursing into `$subkernels`, if any) and replaces every
#' `wrapped_param` named in `axes` whose value has more than one element
#' with its `b`-th element. Parameters not named in `axes` -- and length-1
#' (shared) occurrences of those that are -- are left untouched. `axes` is
#' the contract resolved by [.resolve_axes()]. Used by [batch_kernel()] and
#' [input_specific_param_kernel()].
#'
#' @param kernel A kernel object.
#' @param b A positive integer index.
#' @param axes A character vector of hyperparameter names (the resolved
#'   batch axes).
#' @return A new kernel object, with every axis hyperparameter sliced at
#'   index `b`.
#' @keywords internal
slice_kernel <- function(kernel, b, axes) {
  for (name in names(kernel)) {
    field <- kernel[[name]]
    if (inherits(field, "wrapped_param") && name %in% axes &&
          length(field$wrapped) > 1) {
      kernel[[name]]$wrapped <- field$wrapped[b]
    }
  }
  if (!is.null(kernel$subkernels)) {
    kernel$subkernels <- lapply(kernel$subkernels, slice_kernel, b = b, axes = axes)
  }
  kernel
}

#' @keywords internal
.slice_3d_input <- function(x, b, argument_name = "batch_over_inputs") {
  if (length(dim(x)) != 3) {
    stop(
      "With `", argument_name, " = TRUE`, inputs must be 3-dimensional arrays ",
      "of shape (batch_size, n_points, n_dims).",
      call. = FALSE
    )
  }
  matrix(x[b, , ], nrow = dim(x)[2], ncol = dim(x)[3])
}

#' @keywords internal
.batch_slice_input <- function(x, b, batch_over_inputs) {
  if (!batch_over_inputs) {
    return(x)
  }
  .slice_3d_input(x, b, "batch_over_inputs")
}

#' Add batch handling to a kernel
#'
#' @description
#' Wraps `inner` so that [evaluate()] produces `batch_size` covariance
#' matrices at once, one per batch element, stacked into an
#' `nrow(x1) x nrow(x2) x batch_size` array.
#'
#' Each batch element can use its own hyperparameter values: build `inner`
#' with a hyperparameter of length `batch_size` (one value per batch
#' element) rather than a single scalar, e.g. `se_kernel(length_scale =
#' c(1, 2, 3))` for a batch of 3. A hyperparameter left as a single scalar
#' is shared across every batch element (see [slice_kernel()]).
#'
#' Which hyperparameters vary along the batch axis is an explicit contract
#' (see [.resolve_axes()]): by default (`batch_axes = NULL`) it is
#' auto-detected -- every sliceable multi-valued hyperparameter is assumed
#' to be batched, and must then have exactly `batch_size` values (any other
#' length is an error, not a silent recycling). Name the batched
#' hyperparameters explicitly (`batch_axes = "length_scale"`) whenever
#' auto-detection is ambiguous, most notably for *nested* batching: a
#' hyperparameter not named in `batch_axes` is left untouched for another
#' (inner or outer) `batch_kernel()` to handle, which is what lets each
#' nesting level own its axis independently, with different sizes per
#' level. The resolved axes are validated again at every [evaluate()] call,
#' so a [kupdate()] that changes a value's length cannot silently
#' desynchronise the contract.
#'
#' If `batch_over_inputs = TRUE`, `x1`/`x2` must themselves be 3-dimensional
#' arrays of shape `(batch_size, n_points, n_dims)` -- a different set of
#' input points per batch element, regardless of whether hyperparameters
#' are batched. Defaults to `FALSE` (every batch element sees the same
#' `x1`/`x2`), the more common case for batched hyperparameters.
#'
#' A `batch_kernel()` can itself be wrapped inside another `batch_kernel()`,
#' to handle more than one batch dimension at once; declare each level's
#' `batch_axes` explicitly in that case. This works well when different
#' hyperparameters vary along different axes (e.g. `length_scale` along one
#' level, `noise` along another). When the *same* hyperparameter must vary
#' along the full Cartesian product of two or more axes at once (e.g. a
#' `variance` that needs an independent value for every (cluster, output)
#' combination), nesting cannot express it -- each level's [slice_kernel()]
#' collapses the hyperparameter to a scalar before the inner level ever sees
#' it. Use **multi-axis batching** instead: pass a named integer vector as
#' `batch_size` (e.g. `c(K = 3, O = 2)`) and a named list as `batch_axes`
#' (hyperparameter name -> character vector of the axis names, a subset of
#' `names(batch_size)`, it varies along), all in a single `batch_kernel()`
#' call. Auto-detection is not available in this mode -- `batch_axes` must
#' be given explicitly, since a hyperparameter's length alone cannot
#' disambiguate which axis (or combination of axes) it belongs to. A
#' hyperparameter declared to vary along axes `c(a1, a2, ...)` must have
#' exactly `prod(batch_size[c(a1, a2, ...)])` values (or a single shared
#' one), stored in column-major order over those axes' sizes -- i.e. `a1`
#' is the fastest-varying axis, exactly the order `array(seq_len(prod(sizes)),
#' dim = sizes)` would fill. A hyperparameter can also vary along only a
#' subset of the wrapper's axes (broadcasting across the others), or not be
#' named in `batch_axes` at all (fully shared).
#'
#' Note: the result is an array of shape `(n1, n2, batch_size)` (the batch
#' dimension(s) last, in `names(batch_size)` order for multi-axis) -- the
#' natural shape produced by `simplify2array()`, R's idiomatic way to stack
#' per-element matrices; use `aperm()` to reorder dimensions if a different
#' order is needed. Nested (single-axis) batching appends one more trailing
#' dimension per level: `(n1, n2, inner_batch_size, outer_batch_size)`;
#' multi-axis batching appends one dimension per named axis in a single
#' level: `(n1, n2, batch_size["a1"], batch_size["a2"], ...)`.
#'
#' @param inner A kernel object, or a numeric value (coerced via
#'   [as_kernel()]).
#' @param batch_size A positive integer (single-axis batching: the number
#'   of batch elements), or a named integer vector with 2+ unique,
#'   non-empty names (multi-axis batching: one size per named axis, e.g.
#'   `c(K = 3, O = 2)`).
#' @param batch_over_inputs Logical, whether `x1`/`x2` are themselves
#'   batched 3-dimensional arrays. Defaults to `FALSE`. In multi-axis mode,
#'   the leading dimension of `x1`/`x2` must be `prod(batch_size)`, indexed
#'   in the same column-major order as the output array.
#' @param batch_axes Single-axis mode: `NULL` (auto-detect, the default),
#'   or a character vector naming the hyperparameters that vary along this
#'   batch axis. Multi-axis mode: a required named list (hyperparameter
#'   name -> character vector of axis names); auto-detection is not
#'   supported.
#'
#' @return An object of class `batch_kernel`/`wrapper_kernel`/`kernel`.
#' @export
#' @examples
#' k <- batch_kernel(se_kernel(length_scale = c(1, 2, 3)), batch_size = 3)
#' evaluate(k, c(0, 1, 2)) # array of shape (3, 3, 3)
#'
#' # Nested batching, one axis per level (declared explicitly):
#' k2 <- se_kernel(length_scale = c(.5, 1.5)) + white_noise_kernel(noise = c(0.1, 1, 2, 3))
#' bk <- batch_kernel(k2, batch_size = 2, batch_axes = "length_scale")
#' bbk <- batch_kernel(bk, batch_size = 4, batch_axes = "noise")
#' dim(evaluate(bbk, c(1, 2, 3))) # (3, 3, 2, 4)
#'
#' # Multi-axis batching: `variance` varies over the full 3x2 = 6
#' # (cluster, output) grid in a single wrapper (column-major: cluster
#' # fastest-varying, so values 1:6 map to (K=1,O=1),(K=2,O=1),(K=3,O=1),
#' # (K=1,O=2),(K=2,O=2),(K=3,O=2)):
#' k3 <- variance_kernel(variance = 1:6) * se_kernel(length_scale = 1)
#' mk <- batch_kernel(k3, batch_size = c(K = 3, O = 2),
#'                     batch_axes = list(variance = c("K", "O")))
#' dim(evaluate(mk, c(0, 1, 2))) # (3, 3, 3, 2)
batch_kernel <- function(inner, batch_size, batch_over_inputs = FALSE,
                          batch_axes = NULL) {
  inner <- as_kernel(inner)

  if (length(batch_size) > 1) {
    batch_size <- .validate_multi_batch_size(batch_size, "batch_size")
    batch_axes <- .resolve_axes_multi(inner, batch_size, batch_axes,
                                       "batch_axes", "batch_size")
    kernel <- structure(
      list(
        subkernels = list(inner),
        batch_size = batch_size,
        batch_over_inputs = isTRUE(batch_over_inputs),
        batch_axes = batch_axes
      ),
      class = c("batch_kernel", "wrapper_kernel", "kernel")
    )
    validate_kernel(kernel)
    return(kernel)
  }

  batch_size <- as.integer(batch_size)
  if (length(batch_size) != 1 || is.na(batch_size) || batch_size < 1) {
    stop(
      "`batch_size` must be a single positive integer, or a named integer ",
      "vector of 2+ axis sizes (e.g. `c(K = 3, O = 2)`) for multi-axis ",
      "batching.",
      call. = FALSE
    )
  }
  batch_axes <- .resolve_axes(inner, batch_size, batch_axes,
                               "batch_axes", "batch_size")

  kernel <- structure(
    list(
      subkernels = list(inner),
      batch_size = batch_size,
      batch_over_inputs = isTRUE(batch_over_inputs),
      batch_axes = batch_axes
    ),
    class = c("batch_kernel", "wrapper_kernel", "kernel")
  )
  validate_kernel(kernel)
  kernel
}

#' @export
evaluate.batch_kernel <- function(kernel, x1, x2 = NULL,
                                   na_action = "propagate", mask_variance = 1, ...,
                                   .self_covariance = NULL) {
  inner <- kernel$subkernels[[1]]
  x2_was_null <- is.null(x2)
  if (is.null(.self_covariance)) {
    .self_covariance <- x2_was_null || identical(x1, x2)
  }

  # When batch_over_inputs = FALSE, every slice sees the same x1/x2 -- the
  # distance matrix of any stationary/dot-product leaf is then identical
  # across the K slices (only shape() depends on the per-slice
  # hyperparameter values). Building one plan (R/plan.R) and calling
  # eval_plan() per slice instead of evaluate() skips that recomputation
  # entirely for every cacheable leaf, with an automatic, per-subtree
  # fallback to evaluate() for anything the plan doesn't know how to cache
  # -- see dev/ameliorations/README.md, section C (finding 3a).
  if (!kernel$batch_over_inputs) {
    plan <- precompute(inner, x1, x2)
    evaluate_slice <- function(sliced, x1_b, x2_b) {
      eval_plan(plan, sliced, na_action = na_action, mask_variance = mask_variance, ...,
                .self_covariance = .self_covariance)
    }
  } else {
    evaluate_slice <- function(sliced, x1_b, x2_b) {
      evaluate(sliced, x1_b, x2_b, na_action = na_action, mask_variance = mask_variance, ...,
               .self_covariance = .self_covariance)
    }
  }

  if (length(kernel$batch_size) > 1) {
    # Re-validate the axis contract: a kupdate() since construction may
    # have changed a hyperparameter's length.
    axis_spec <- .resolve_axes_multi(inner, kernel$batch_size, kernel$batch_axes,
                                      "batch_axes", "batch_size")
    n_total <- prod(kernel$batch_size)

    results <- lapply(seq_len(n_total), function(idx) {
      combo <- stats::setNames(
        as.integer(arrayInd(idx, .dim = as.integer(kernel$batch_size))[1, ]),
        names(kernel$batch_size)
      )
      sliced <- slice_kernel_nd(inner, combo, axis_spec, kernel$batch_size)
      x1_b <- .batch_slice_input(x1, idx, kernel$batch_over_inputs)
      x2_b <- if (x2_was_null) NULL else .batch_slice_input(x2, idx, kernel$batch_over_inputs)
      evaluate_slice(sliced, x1_b, x2_b)
    })

    # Built directly via array()/unlist() rather than simplify2array()+dim<-:
    # simplify2array() drops the dim attribute entirely (falling back to a
    # bare vector) whenever every per-combo result has length 1 -- the
    # single-input-point self-covariance case -- which would silently
    # reshape the flat vector into just (batch_size...) instead of
    # (1, 1, batch_size...), losing the covariance-matrix shape. Deriving
    # n1/n2 from the first result via as.matrix() (rather than from the
    # combined array's own, unreliable post-hoc dim) sidesteps that case.
    first <- as.matrix(results[[1]])
    return(array(unlist(results), dim = c(nrow(first), ncol(first), unname(kernel$batch_size))))
  }

  # Re-validate the axis contract: a kupdate() since construction may have
  # changed a hyperparameter's length.
  axes <- .resolve_axes(inner, kernel$batch_size, kernel$batch_axes,
                         "batch_axes", "batch_size")

  results <- lapply(seq_len(kernel$batch_size), function(b) {
    sliced <- slice_kernel(inner, b, axes)
    x1_b <- .batch_slice_input(x1, b, kernel$batch_over_inputs)
    x2_b <- if (x2_was_null) NULL else .batch_slice_input(x2, b, kernel$batch_over_inputs)
    evaluate_slice(sliced, x1_b, x2_b)
  })

  # simplify2array() builds the right shape whenever each per-batch result
  # has length > 1 -- (n1, n2, batch_size) for a matrix (dense engine),
  # (n, batch_size) for a vector (diagonal engine) -- verified empirically,
  # not just assumed. It only misbehaves for the single-input-point
  # self-covariance case, where every result has length exactly 1 (a 1x1
  # matrix, or a length-1 diagonal vector): it then drops the dim attribute
  # entirely, silently returning a bare length-batch_size vector instead of
  # a (1, 1, batch_size) array (dense) or (1, batch_size) matrix (diagonal).
  # That broke callers expecting an array/matrix unconditionally -- e.g.
  # block_diag_kernel() with nb_blocks = 1 on a single point, which errored
  # outright trying to derive nrow()/ncol() from a bare vector's NULL dim.
  # Handled explicitly here, only for that degenerate length-1 case.
  first <- results[[1]]
  if (length(first) > 1) {
    return(simplify2array(results))
  }
  if (is.matrix(first)) {
    array(unlist(results), dim = c(1L, 1L, kernel$batch_size))
  } else {
    matrix(unlist(results), nrow = 1L, ncol = kernel$batch_size)
  }
}

#' @export
format.batch_kernel <- function(x, ...) {
  # Just the inner kernel's format(): the batch information lives in the
  # shape of its hyperparameter values, not in a separate notation.
  format(x$subkernels[[1]])
}

#' Per-input hyperparameters
#'
#' @description
#' Wraps `inner` so that every point of `x1` gets its own hyperparameter
#' value, rather than a single value shared across all of `x1`. Build
#' `inner` with a hyperparameter of length `input_size` to vary it per
#' point of `x1` (see [slice_kernel()]); a value left as a single scalar is
#' shared across all points of `x1`. A typical use is heteroscedastic noise:
#' `input_specific_param_kernel(white_noise_kernel(noise = c(0.1, 0.5, 0.2)),
#' input_size = 3)` gives each of the 3 points of `x1` its own noise level.
#'
#' Parameters are only varied along `x1`; `x2` is always evaluated as a
#' whole against each sliced point of `x1` in turn, so take care when `x1`
#' and `x2` differ.
#'
#' **Symmetry caveat**: in the auto-covariance case (`x2 = NULL`), entry
#' `(i, j)` is computed with point `i`'s hyperparameter value and entry
#' `(j, i)` with point `j`'s -- the result is only a symmetric (and
#' positive semi-definite) covariance matrix if the varied hyperparameter
#' acts on the *diagonal only*, as [white_noise_kernel()]'s `noise` does.
#' Varying a cross-covariance hyperparameter (e.g. a `length_scale`)
#' per point produces an asymmetric matrix, which is not a valid GP
#' covariance.
#'
#' Which hyperparameters vary per input point is the same explicit contract
#' as [batch_kernel()]'s (see [.resolve_axes()]): auto-detected by default
#' (`param_axes = NULL`, every sliceable multi-valued hyperparameter must
#' then have exactly `input_size` values), or named explicitly via
#' `param_axes`, leaving other multi-valued hyperparameters to other
#' wrappers' axes.
#'
#' @param inner A kernel object, or a numeric value (coerced via
#'   [as_kernel()]).
#' @param input_size A positive integer: the expected number of points in
#'   `x1` (an error is raised if `nrow(x1)` does not match).
#' @param param_axes `NULL` (auto-detect, the default), or a character
#'   vector naming the hyperparameters that vary per input point.
#'
#' @return An object of class `input_specific_param_kernel`/
#'   `wrapper_kernel`/`kernel`.
#' @export
#' @examples
#' k <- input_specific_param_kernel(
#'   white_noise_kernel(noise = c(0.1, 0.5, 0.2)),
#'   input_size = 3
#' )
#' evaluate(k, c(0, 1, 2))
input_specific_param_kernel <- function(inner, input_size, param_axes = NULL) {
  input_size <- as.integer(input_size)
  if (length(input_size) != 1 || input_size < 1) {
    stop("`input_size` must be a single positive integer.", call. = FALSE)
  }
  inner <- as_kernel(inner)
  param_axes <- .resolve_axes(inner, input_size, param_axes,
                               "param_axes", "input_size")

  kernel <- structure(
    list(
      subkernels = list(inner),
      input_size = input_size,
      param_axes = param_axes
    ),
    class = c("input_specific_param_kernel", "wrapper_kernel", "kernel")
  )
  validate_kernel(kernel)
  kernel
}

#' @export
evaluate.input_specific_param_kernel <- function(kernel, x1, x2 = NULL,
                                                  na_action = c("propagate", "mask"),
                                                  mask_variance = 1, ...,
                                                  .self_covariance = NULL) {
  na_action <- match.arg(na_action)
  inner <- kernel$subkernels[[1]]
  X1 <- as_matrix_input(x1)
  if (nrow(X1) != kernel$input_size) {
    stop(
      "nrow(x1) (", nrow(X1), ") does not match `input_size` (", kernel$input_size, ").",
      call. = FALSE
    )
  }
  # The self-covariance question can only be answered *here*, before x2 is
  # materialised: each per-row sub-call below compares a 1-row x1 against
  # the full x2 and would always (wrongly) conclude "not a self-covariance".
  if (is.null(.self_covariance)) {
    .self_covariance <- is.null(x2) || identical(x1, x2)
  }
  if (is.null(x2)) {
    x2 <- x1
  }
  axes <- .resolve_axes(inner, kernel$input_size, kernel$param_axes,
                         "param_axes", "input_size")

  # Sub-calls always propagate NA: masking needs the full-matrix geometry
  # (the (i, i) diagonal position of a masked point), which a single 1 x n2
  # row does not have -- it is applied once, below, on the assembled matrix.
  rows <- lapply(seq_len(kernel$input_size), function(i) {
    sliced <- slice_kernel(inner, i, axes)
    evaluate(sliced, X1[i, , drop = FALSE], x2, na_action = "propagate", ...,
             .self_covariance = FALSE)
  })
  result <- do.call(rbind, rows)

  if (na_action == "mask") {
    X2 <- as_matrix_input(x2)
    na1 <- apply(X1, 1, anyNA)
    na2 <- apply(X2, 1, anyNA)
    result <- apply_na_mask(result, na1, na2, .self_covariance, mask_variance)
  }
  result
}

#' @export
format.input_specific_param_kernel <- function(x, ...) {
  format(x$subkernels[[1]])
}

# block_diag_kernel()/block_kernel(): build block covariance matrices.
# block_diag_kernel() is the easy case -- exactly batch_kernel() plus
# assembling the per-block matrices block-diagonally, since each block is
# independent (no interaction term needed between block i and block j).
# block_kernel() is harder: it also computes the *cross*-block terms, which
# only make sense if `inner` can read a *pair* of hyperparameter values
# (one per block side) -- see pair_slice_kernel() below.

#' Block-diagonal covariance matrix
#'
#' @description
#' Wraps `inner` so that [evaluate()] produces a block-diagonal matrix of
#' `nb_blocks` independent blocks (off-diagonal blocks are `0`) -- exactly
#' [batch_kernel()], with the resulting `(n1, n2, nb_blocks)` array
#' reassembled into a single `(n1 * nb_blocks) x (n2 * nb_blocks)` matrix.
#' Each block may have its own hyperparameter values the same way
#' [batch_kernel()]'s do (see [slice_kernel()]).
#'
#' @param inner A kernel object, or a numeric value (coerced via
#'   [as_kernel()]).
#' @param nb_blocks A positive integer: the number of blocks.
#' @param batch_over_inputs Logical, whether `x1`/`x2` are themselves
#'   per-block 3-dimensional arrays. Defaults to `FALSE`. See
#'   [batch_kernel()].
#' @param batch_axes `NULL` (auto-detect, the default), or a character
#'   vector naming the hyperparameters that vary per block. See
#'   [batch_kernel()].
#'
#' @return An object of class `block_diag_kernel`/`wrapper_kernel`/`kernel`.
#' @export
#' @examples
#' k <- block_diag_kernel(se_kernel(length_scale = c(1, 2)), nb_blocks = 2)
#' evaluate(k, c(0, 1)) # 4x4, block-diagonal
block_diag_kernel <- function(inner, nb_blocks, batch_over_inputs = FALSE,
                               batch_axes = NULL) {
  nb_blocks_check <- as.integer(nb_blocks)
  if (length(nb_blocks_check) != 1 || is.na(nb_blocks_check) || nb_blocks_check < 1) {
    stop(
      "`nb_blocks` must be a single positive integer -- block_diag_kernel() ",
      "has no multi-axis mode (unlike batch_kernel()), since ",
      "evaluate.block_diag_kernel() only knows how to reassemble a single ",
      "(n1, n2, nb_blocks) array into a block-diagonal matrix.",
      call. = FALSE
    )
  }
  kernel <- structure(
    list(
      subkernels = list(batch_kernel(inner, batch_size = nb_blocks,
                                      batch_over_inputs = batch_over_inputs,
                                      batch_axes = batch_axes)),
      nb_blocks = as.integer(nb_blocks)
    ),
    class = c("block_diag_kernel", "wrapper_kernel", "kernel")
  )
  validate_kernel(kernel)
  kernel
}

#' @export
evaluate.block_diag_kernel <- function(kernel, x1, x2 = NULL, ...) {
  blocks <- evaluate(kernel$subkernels[[1]], x1, x2, ...) # (n1, n2, nb_blocks)
  n1 <- dim(blocks)[1]
  n2 <- dim(blocks)[2]
  nb_blocks <- dim(blocks)[3]

  result <- matrix(0, nrow = n1 * nb_blocks, ncol = n2 * nb_blocks)
  for (b in seq_len(nb_blocks)) {
    rows <- ((b - 1) * n1 + 1):(b * n1)
    cols <- ((b - 1) * n2 + 1):(b * n2)
    result[rows, cols] <- blocks[, , b]
  }
  result
}

#' @export
format.block_diag_kernel <- function(x, ...) {
  paste0("BlockDiag", format(x$subkernels[[1]]))
}

# block_kernel()'s cross-block terms require `inner` to read a *pair* of
# hyperparameter values (one per block side) rather than a single shared
# scalar -- found by re-examining feature_kernel(), which already does
# exactly this for its length_scale/variance pairs (kern_length_scale[1]/
# [2]). pair_slice_kernel() generalises slice_kernel() to extract such a
# pair (i, j) instead of a single index b.

#' Extract a pair of batch/block elements of a kernel's hyperparameters
#'
#' @description
#' Like [slice_kernel()], but replaces every `wrapped_param` named in
#' `axes` (with more than one value) with the **pair** of its `i`-th and
#' `j`-th elements (a length-2 vector), rather than a single element. Used
#' by [block_kernel()] to build the cross-block (block `i`, block `j`)
#' term, which needs both blocks' hyperparameter values available at once
#' (e.g. [feature_kernel()]'s paired `length_scale`/`variance`).
#'
#' @param kernel A kernel object.
#' @param i,j Positive integer indices.
#' @param axes A character vector of hyperparameter names (the resolved
#'   block axes).
#' @return A new kernel object, with every axis hyperparameter replaced by
#'   its `c(i, j)` pair.
#' @keywords internal
pair_slice_kernel <- function(kernel, i, j, axes) {
  for (name in names(kernel)) {
    field <- kernel[[name]]
    if (inherits(field, "wrapped_param") && name %in% axes &&
          length(field$wrapped) > 1) {
      kernel[[name]]$wrapped <- field$wrapped[c(i, j)]
    }
  }
  if (!is.null(kernel$subkernels)) {
    kernel$subkernels <- lapply(kernel$subkernels, pair_slice_kernel,
                                 i = i, j = j, axes = axes)
  }
  kernel
}

#' Block covariance matrix, with cross-block terms
#'
#' @description
#' Wraps `inner` so that [evaluate()] produces a full `(n1 * nb_blocks) x
#' (n2 * nb_blocks)` block matrix, including the cross-block (block `i`,
#' block `j`, `i != j`) terms -- unlike [block_diag_kernel()], which leaves
#' them at `0`.
#'
#' Computing a cross-block term requires `inner` to read a **pair** of
#' hyperparameter values at once (one per block side), via
#' [pair_slice_kernel()]: build `inner` with a hyperparameter of length
#' `nb_blocks` to vary it per block, exactly as for [batch_kernel()] (see
#' [slice_kernel()]).
#'
#' Most kernels are not designed to interpret a *pair* of values per
#' hyperparameter ([se_kernel()], for instance, expects a single scalar
#' `length_scale`) -- using one of those with per-block hyperparameter
#' values would silently produce a mathematically meaningless result.
#' [feature_kernel()] is, deliberately, designed for exactly this (it
#' already reads a length-2 `length_scale`/`variance` pair, one value per
#' "side"). To catch misuse early rather than silently computing garbage,
#' `block_kernel()` raises an error if `inner` has a per-block
#' hyperparameter (length `nb_blocks`) but is not explicitly marked
#' block-compatible (`attr(inner, "block_compatible")` -- see
#' [feature_kernel()]'s constructor). A kernel with only shared (length-1)
#' hyperparameters needs no such marker: every block is then identical, no
#' pairing is ever attempted.
#'
#' Only the upper triangle of blocks is actually computed; the rest is
#' filled by symmetry (block `(j, i)` is the transpose of block `(i, j)`),
#' since `inner`'s formula is symmetric under swapping which side is
#' "left"/"right" together with `x1`/`x2` -- this algebraic shortcut is
#' expressed here as a dedicated mapping primitive, `map_block_symmetric()`
#' (layer 2 of the architecture), rather than inline in `evaluate()`.
#'
#' @param inner A kernel object, or a numeric value (coerced via
#'   [as_kernel()]).
#' @param nb_blocks A positive integer: the number of blocks.
#' @param batch_over_inputs Logical, whether `x1`/`x2` are themselves
#'   per-block 3-dimensional arrays. Defaults to `FALSE`. See
#'   [batch_kernel()].
#' @param block_axes `NULL` (auto-detect, the default), or a character
#'   vector naming the hyperparameters that vary per block. Same contract
#'   as [batch_kernel()]'s `batch_axes` (see [.resolve_axes()]).
#'
#' @return An object of class `block_kernel`/`wrapper_kernel`/`kernel`.
#' @export
#' @examples
#' k <- block_kernel(
#'   feature_kernel(length_scale = c(1, 2), length_scale_u = 1, variance = c(1, 1)),
#'   nb_blocks = 2
#' )
#' evaluate(k, c(0, 1)) # 4x4, with non-zero cross-block terms
block_kernel <- function(inner, nb_blocks, batch_over_inputs = FALSE,
                          block_axes = NULL) {
  nb_blocks <- as.integer(nb_blocks)
  if (length(nb_blocks) != 1 || nb_blocks < 1) {
    stop("`nb_blocks` must be a single positive integer.", call. = FALSE)
  }
  inner <- as_kernel(inner)
  block_axes <- .resolve_axes(inner, nb_blocks, block_axes,
                               "block_axes", "nb_blocks")
  if (length(block_axes) > 0 && !isTRUE(attr(inner, "block_compatible"))) {
    stop(
      "`inner` must be explicitly marked `block_compatible` (e.g. feature_kernel()) ",
      "to be used with block_kernel() when its hyperparameters vary per block. ",
      "Standard kernels (e.g. se_kernel()) do not know how to interpret a pair of ",
      "values for a single hyperparameter.",
      call. = FALSE
    )
  }

  kernel <- structure(
    list(
      subkernels = list(inner),
      nb_blocks = nb_blocks,
      batch_over_inputs = isTRUE(batch_over_inputs),
      block_axes = block_axes
    ),
    class = c("block_kernel", "wrapper_kernel", "kernel")
  )
  validate_kernel(kernel)
  kernel
}

#' Build a symmetric block matrix, computing only the upper triangle of blocks
#'
#' @description
#' The mapping primitive backing [block_kernel()]: computes block `(i, j)`
#' for every `i <= j` via [pair_slice_kernel()] + [evaluate()], and fills
#' block `(j, i)` as its transpose -- exploiting the symmetry of `kernel`'s
#' formula under swapping sides.
#'
#' @param kernel A kernel object (the `inner` of a `block_kernel`).
#' @param X1,X2 Numeric matrices, or 3-dimensional arrays if
#'   `block_over_inputs = TRUE`.
#' @param nb_blocks A positive integer.
#' @param block_over_inputs Logical.
#' @param ... Passed on to [evaluate()].
#' @return A numeric matrix of dimension `(nrow(X1) * nb_blocks) x (nrow(X2) * nb_blocks)`
#'   (or, with `block_over_inputs = TRUE`, `(dim(X1)[2] * nb_blocks) x (dim(X2)[2] * nb_blocks)`).
#' @keywords internal
map_block_symmetric <- function(kernel, X1, X2, nb_blocks, block_over_inputs,
                                 axes, self_covariance, ...) {
  block_grid <- vector("list", nb_blocks)
  for (i in seq_len(nb_blocks)) {
    block_grid[[i]] <- vector("list", nb_blocks)
  }

  for (i in seq_len(nb_blocks)) {
    for (j in i:nb_blocks) {
      sliced <- pair_slice_kernel(kernel, i, j, axes)
      x1_b <- if (block_over_inputs) .slice_3d_input(X1, i, "batch_over_inputs") else X1
      x2_b <- if (block_over_inputs) .slice_3d_input(X2, j, "batch_over_inputs") else X2
      # Only the diagonal blocks (i == j) of a self-covariance are
      # themselves self-covariances; a cross block (i != j) covaries two
      # *different* processes over the same points -- its diagonal must not
      # receive the self-variance treatment (e.g. na_action = "mask"'s
      # mask_variance).
      block <- evaluate(sliced, x1_b, x2_b, ...,
                        .self_covariance = self_covariance && i == j)

      block_grid[[i]][[j]] <- block
      if (j != i) {
        block_grid[[j]][[i]] <- t(block)
      }
    }
  }

  do.call(rbind, lapply(seq_len(nb_blocks), function(i) do.call(cbind, block_grid[[i]])))
}

#' @export
evaluate.block_kernel <- function(kernel, x1, x2 = NULL, ...,
                                   .self_covariance = NULL) {
  if (is.null(.self_covariance)) {
    .self_covariance <- is.null(x2) || identical(x1, x2)
  }
  if (is.null(x2)) {
    x2 <- x1
  }
  # Re-validate the axis contract (see evaluate.batch_kernel).
  axes <- .resolve_axes(kernel$subkernels[[1]], kernel$nb_blocks,
                         kernel$block_axes, "block_axes", "nb_blocks")
  map_block_symmetric(kernel$subkernels[[1]], x1, x2, kernel$nb_blocks,
                       kernel$batch_over_inputs, axes, .self_covariance, ...)
}

#' @export
format.block_kernel <- function(x, ...) {
  paste0("Block", format(x$subkernels[[1]]))
}
