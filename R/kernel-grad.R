# Analytic gradients mirror the evaluate() pipeline one layer at a time
# (see dev/ameliorations/README.md, backlog item 2b): shape_grad() sits next
# to shape() (R/pairwise.R's fast vectorised path), pairwise_matrix_grad()/
# pairwise_diag_grad() sit next to pairwise_matrix()/pairwise_diag()
# (R/distance-matrix.R), engine_grad() sits next to engine_eval()
# (R/engines.R), and kernel_grad() sits next to evaluate() (R/evaluate.R).
# kernel_grad_raw() is the recursive tree-walk underneath kernel_grad(),
# analogous to flatten_params() underneath get_params() -- it returns one
# entry per trainable hyperparameter, in the same own-fields-then-subkernels
# order flatten_params()/.flatten_params_with_class() use, keyed by bare
# field name (duplicates allowed); kernel_grad() disambiguates once at the
# end, the same way get_params()/get_trainable_params() do.

#' @keywords internal
.kernel_grad_unsupported <- function(kernel) {
  stop(
    "kernel_grad() is not implemented for `", class(kernel)[1], "` -- analytic ",
    "gradients are only available for the base kernels (se_kernel, ",
    "matern12_kernel/matern32_kernel/matern52_kernel, periodic_kernel, ",
    "rational_quadratic_kernel, white_noise_kernel, feature_kernel, ",
    "linear_kernel, affine_kernel, polynomial_kernel, sigmoid_kernel, ",
    "constant_kernel, variance_kernel), their sum_kernel()/product_kernel() ",
    "compositions, exp_kernel(), log_kernel(), active_dims_kernel(), and ",
    "ard_kernel() (see ?kernel_grad for ard_kernel()'s extra requirements). ",
    "Not supported: batch_kernel(), block_kernel(), block_diag_kernel(), ",
    "input_specific_param_kernel().",
    call. = FALSE
  )
}

#' Elementwise gradient of a stationary/dot-product kernel's `shape()`
#'
#' @description
#' The derivative counterpart of [shape()]: for each of `kernel`'s own
#' hyperparameters, the partial derivative of `shape(kernel, d)` with respect
#' to that hyperparameter (in its natural, constrained space), evaluated at
#' `d`. Like `shape()`, this is a purely elementwise scalar transform, so it
#' vectorises over `d` being a value, vector, or matrix for free.
#'
#' A `shape_grad()` method is what lets [pairwise_matrix_grad()]/
#' [pairwise_diag_grad()] give a `stationary_kernel`/`dot_product_kernel` its
#' fast, vectorised gradient path -- the same role `shape()` plays for
#' [pairwise_matrix()]/[pairwise_diag()]. A kernel with no `shape_grad()`
#' method is simply not covered by [kernel_grad()] (see its documentation for
#' the exact list of kernels that are).
#'
#' @param kernel A kernel object.
#' @param d A numeric value, vector, or matrix of distances (as passed to
#'   [shape()]).
#' @return A named list, one entry per own hyperparameter of `kernel`
#'   (matching its field names), each the same shape as `d`.
#' @export
#' @examples
#' k <- se_kernel(length_scale = 2)
#' shape_grad(k, c(0, 1, 4))
shape_grad <- function(kernel, d) {
  UseMethod("shape_grad")
}

#' Elementwise derivative of `shape()` with respect to the distance itself
#'
#' @description
#' Unlike [shape_grad()] (the derivative with respect to one of `kernel`'s
#' own hyperparameters), `shape_grad_d()` differentiates `shape(kernel, d)`
#' with respect to `d` itself. This is the piece [ard_kernel()] needs: it
#' feeds its inner kernel a *transformed* distance (computed on `x /
#' length_scales`), so getting from a `length_scales` gradient back to a
#' hyperparameter gradient means chain-ruling through `shape()` at the
#' distance level, not the hyperparameter level.
#'
#' @param kernel A kernel object.
#' @param d A numeric value, vector, or matrix of distances (as passed to
#'   [shape()]).
#' @return The same shape as `d`.
#' @export
#' @examples
#' k <- se_kernel(length_scale = 2)
#' shape_grad_d(k, c(0, 1, 4))
shape_grad_d <- function(kernel, d) {
  UseMethod("shape_grad_d")
}

#' Gradient of the vectorised pairwise covariance matrix/diagonal
#'
#' @description
#' The derivative counterparts of [pairwise_matrix()]/[pairwise_diag()]:
#' for each of `kernel`'s own hyperparameters, the partial derivative of the
#' covariance matrix (or, for `pairwise_diag_grad()`, its diagonal) with
#' respect to that hyperparameter. `stationary_kernel`/`dot_product_kernel`
#' subclasses get a fast method via [distance_matrix()]/[diag_distance()] +
#' [shape_grad()], the same pattern [pairwise_matrix()]/[pairwise_diag()]
#' use with [shape()]. A handful of kernels without a `shape()` method
#' (`constant_kernel()`, `variance_kernel()`, `feature_kernel()`) get a
#' direct method instead.
#'
#' @param kernel A kernel object.
#' @param X1,X2 Numeric matrices.
#' @return A named list, one entry per own hyperparameter of `kernel`. A
#'   scalar-valued hyperparameter's entry is a numeric matrix (dimension
#'   `nrow(X1)` x `nrow(X2)`, for `pairwise_matrix_grad()`) or vector (length
#'   `nrow(X1)`, for `pairwise_diag_grad()`); a vector-valued hyperparameter
#'   (currently only `feature_kernel()`'s `length_scale`/`variance`) gets a
#'   *list* of one such matrix/vector per component instead, in the same
#'   order as the hyperparameter's own value.
#' @name pairwise_grad
NULL

#' @rdname pairwise_grad
#' @keywords internal
pairwise_matrix_grad <- function(kernel, X1, X2) {
  UseMethod("pairwise_matrix_grad")
}

#' @keywords internal
#' @exportS3Method
pairwise_matrix_grad.default <- function(kernel, X1, X2) {
  .kernel_grad_unsupported(kernel)
}

#' @keywords internal
#' @exportS3Method
pairwise_matrix_grad.stationary_kernel <- function(kernel, X1, X2) {
  if (!has_s3_method("shape_grad", kernel)) {
    return(pairwise_matrix_grad.default(kernel, X1, X2))
  }
  d <- distance_matrix(kernel$distance_function, X1, X2)
  shape_grad(kernel, d)
}

#' @keywords internal
#' @exportS3Method
pairwise_matrix_grad.dot_product_kernel <- function(kernel, X1, X2) {
  if (!has_s3_method("shape_grad", kernel)) {
    return(pairwise_matrix_grad.default(kernel, X1, X2))
  }
  dp <- distance_matrix(kernel$distance_function, X1, X2)
  shape_grad(kernel, dp)
}

#' @rdname pairwise_grad
#' @keywords internal
pairwise_diag_grad <- function(kernel, X1, X2) {
  UseMethod("pairwise_diag_grad")
}

#' @keywords internal
#' @exportS3Method
pairwise_diag_grad.default <- function(kernel, X1, X2) {
  .kernel_grad_unsupported(kernel)
}

#' @keywords internal
#' @exportS3Method
pairwise_diag_grad.stationary_kernel <- function(kernel, X1, X2) {
  if (!has_s3_method("shape_grad", kernel)) {
    return(pairwise_diag_grad.default(kernel, X1, X2))
  }
  d <- diag_distance(kernel$distance_function, X1, X2)
  shape_grad(kernel, d)
}

#' @keywords internal
#' @exportS3Method
pairwise_diag_grad.dot_product_kernel <- function(kernel, X1, X2) {
  if (!has_s3_method("shape_grad", kernel)) {
    return(pairwise_diag_grad.default(kernel, X1, X2))
  }
  dp <- diag_distance(kernel$distance_function, X1, X2)
  shape_grad(kernel, dp)
}

#' Gradient counterpart of an engine
#'
#' @description
#' The derivative counterpart of [engine_eval()]: routes to
#' [pairwise_matrix_grad()] or [pairwise_diag_grad()] depending on `engine`,
#' the same way `engine_eval()` routes to [pairwise_matrix()]/
#' [pairwise_diag()].
#'
#' @param engine An engine object.
#' @param kernel A kernel object.
#' @param X1,X2 Numeric matrices.
#' @return See [pairwise_matrix_grad()]/[pairwise_diag_grad()].
#' @name engine_grad
NULL

#' @rdname engine_grad
#' @export
engine_grad <- function(engine, kernel, X1, X2) {
  UseMethod("engine_grad")
}

#' @export
engine_grad.dense_engine <- function(engine, kernel, X1, X2) {
  pairwise_matrix_grad(kernel, X1, X2)
}

#' @export
engine_grad.diagonal_engine <- function(engine, kernel, X1, X2) {
  if (nrow(X1) != nrow(X2)) {
    stop(
      "diagonal_engine() requires `x1` and `x2` to have the same number of points ",
      "(got ", nrow(X1), " and ", nrow(X2), ").",
      call. = FALSE
    )
  }
  pairwise_diag_grad(kernel, X1, X2)
}

#' Elementwise product across gradient entries that may be plain arrays or
#' lists of arrays (the representation used for a vector-valued
#' hyperparameter, see [pairwise_matrix_grad()]).
#'
#' @keywords internal
.grad_times <- function(a, b) {
  if (is.list(a)) {
    return(lapply(a, .grad_times, b = b))
  }
  if (is.list(b)) {
    return(lapply(b, function(bb) .grad_times(a, bb)))
  }
  a * b
}

#' Recursive tree-walk underneath [kernel_grad()]
#'
#' @description
#' Like [flatten_params()] is to [get_params()], `kernel_grad_raw()` is the
#' undisambiguated recursion underneath [kernel_grad()]: one entry per
#' trainable hyperparameter, keyed by bare field name (duplicates allowed
#' when the same name occurs more than once in the tree), in the same
#' own-fields-then-subkernels order [.flatten_params_with_class()] walks.
#' [kernel_grad()] disambiguates the result once, at the very end -- doing
#' it here instead would disambiguate each subtree in isolation, missing
#' collisions with hyperparameters elsewhere in the tree.
#'
#' @param kernel A kernel object.
#' @param x1,x2 As for [evaluate()] (`x2 = NULL` means the auto-covariance
#'   case).
#' @return A named list (see [pairwise_matrix_grad()] for the shape of each
#'   entry).
#' @keywords internal
kernel_grad_raw <- function(kernel, x1, x2) {
  UseMethod("kernel_grad_raw")
}

#' @keywords internal
#' @exportS3Method
kernel_grad_raw.default <- function(kernel, x1, x2) {
  if (is.null(kernel$engine)) {
    .kernel_grad_unsupported(kernel)
  }
  X1 <- as_matrix_input(x1)
  X2 <- as_matrix_input(if (is.null(x2)) x1 else x2)
  own <- Filter(function(field) inherits(field, "wrapped_param"), kernel)
  trainable_names <- names(Filter(function(p) !isTRUE(p$frozen), own))
  if (length(trainable_names) == 0) {
    return(list())
  }
  engine_grad(kernel$engine, kernel, X1, X2)[trainable_names]
}

#' Analytic gradient of a kernel's covariance matrix
#'
#' @description
#' For each non-frozen hyperparameter of `kernel` (in its natural,
#' constrained space -- the same space [get_trainable_params()] reports),
#' the partial derivative of `evaluate(kernel, x1, x2)` with respect to that
#' hyperparameter. Handing this to an optimiser (e.g. `optim(method =
#' "L-BFGS-B", gr = ...)`) alongside the usual objective avoids the `2p + 1`
#' finite-difference evaluations per step a gradient-free optimiser needs
#' (see `dev/ameliorations/README.md`, backlog item 2b) -- one analytic pass
#' instead.
#'
#' `kernel_grad()` does not compute a likelihood gradient itself: like the
#' rest of keRnel, it stays a general-purpose kernel package (see
#' `dev/ameliorations/README.md`, section D) and leaves the reduction (e.g.
#' the Gaussian-process marginal-likelihood trace formula) to the caller,
#' who combines each returned matrix with `K^{-1}` and the residuals however
#' their model requires.
#'
#' Supported on the base kernels -- [se_kernel()], [matern12_kernel()]/
#' [matern32_kernel()]/[matern52_kernel()], [periodic_kernel()],
#' [rational_quadratic_kernel()], [white_noise_kernel()], [feature_kernel()],
#' [linear_kernel()], [affine_kernel()], [polynomial_kernel()],
#' [sigmoid_kernel()], [constant_kernel()], [variance_kernel()] -- any
#' [sum_kernel()]/[product_kernel()] composition of them (including the
#' compositions [kernel_interactions()] builds); [exp_kernel()] and
#' [log_kernel()]; [active_dims_kernel()]; and [ard_kernel()], provided its
#' wrapped kernel has a [shape()] method and one of `squared_euclidean_distance()`
#' /`euclidean_distance()`/`manhattan_distance()`/`dot_product()` as its
#' `distance_function` (a clear error names the unsupported kernel/distance
#' otherwise).
#'
#' Not supported on the remaining wrapper kernels -- [batch_kernel()],
#' [block_kernel()], [block_diag_kernel()], [input_specific_param_kernel()]
#' -- which slice hyperparameters per batch/block/input element rather than
#' compute a single `shape()`/`pairwise_matrix()` formula. Calling
#' `kernel_grad()` on one of these, or on a composition that contains one,
#' raises a clear error naming the unsupported class.
#'
#' Unlike [evaluate()], `kernel_grad()` has no `na_action`: `x1`/`x2` must
#' be free of missing/non-finite values.
#'
#' @param kernel A kernel object (see the "Supported kernels" list above).
#' @param x1 A numeric vector or matrix of input points (see
#'   [as_matrix_input()]).
#' @param x2 An optional second set of input points; defaults to `x1`.
#' @return A named list, one entry per non-frozen hyperparameter -- same
#'   names as [get_trainable_params()] (disambiguated the same way when a
#'   bare name collides elsewhere in the tree). Each entry has the shape
#'   described in [pairwise_matrix_grad()].
#' @export
#' @examples
#' k <- se_kernel(length_scale = 2)
#' x <- c(0, 1, 3)
#' g <- kernel_grad(k, x)
#' g$length_scale
#'
#' # A sum_kernel()/product_kernel() composition works the same way:
#' k2 <- se_kernel(length_scale = 2) + white_noise_kernel(noise = 0.1)
#' names(kernel_grad(k2, x))
kernel_grad <- function(kernel, x1, x2 = NULL) {
  # Reuses evaluate()'s own validation (shape/engine mismatches across a
  # composite tree, .check_scalar_params()...) rather than duplicating it --
  # the result itself is discarded.
  evaluate(kernel, x1, x2)

  raw <- kernel_grad_raw(kernel, x1, x2)

  fp <- .flatten_params_with_class(kernel)
  trainable_mask <- !vapply(fp$params, function(p) isTRUE(p$frozen), logical(1))
  # kernel_grad_raw() walks the same own-fields-then-subkernels order as
  # .flatten_params_with_class(), so the two stay aligned entry-for-entry.
  names(raw) <- .disambiguate_param_names(names(raw), fp$classes[trainable_mask])
  raw
}

#' Analytic gradient of a kernel's covariance matrix, in free (unconstrained) space
#'
#' @description
#' Like [kernel_grad()], but chain-ruled through each hyperparameter's
#' [parametrisation] into the *free* (wrapped, unconstrained) space
#' [get_free_params()]/[set_free_params()] work in, via [unwrap_grad()]
#' (`d(unwrap(p, w))/dw`).
#'
#' The result is unnamed and strictly positional, in exactly the tree-walk
#' order [get_free_params()] uses -- including splitting a vector-valued
#' hyperparameter into one entry per component, the same way
#' [get_free_params()]'s `unlist()` does. This is what makes it safe to
#' zip element-by-element against [get_free_params()]'s own output, or to
#' hand a reduction of it straight to [set_free_params()]'s coordinate
#' system, *even when two sub-kernels share a bare hyperparameter name* --
#' [kernel_grad()]/[get_trainable_params()]'s names disambiguate that case
#' for display, but [kupdate()] cannot address the two occurrences
#' independently by name; working positionally in free space sidesteps the
#' issue entirely, the same way [get_free_params()]/[set_free_params()]
#' already do for the *values* themselves.
#'
#' A worked example: an `optim(method = "L-BFGS-B")` M-step across a
#' Gaussian-process log-likelihood, without box constraints (positivity is
#' already guaranteed by each hyperparameter's parametrisation) --
#'
#' ```r
#' gr <- function(theta, kernel, X, alpha, K_inv) {
#'   k <- set_free_params(kernel, theta)
#'   grads <- kernel_grad_free(k, X)
#'   vapply(grads, function(dK) {
#'     -0.5 * sum(diag((alpha %*% t(alpha) - K_inv) %*% dK))
#'   }, numeric(1))
#' }
#' optim(get_free_params(kernel), fn, gr = gr, method = "L-BFGS-B", kernel = kernel, X = X, ...)
#' ```
#'
#' @param kernel A kernel object (same support as [kernel_grad()]).
#' @param x1,x2 As for [kernel_grad()].
#' @return An unnamed list, `length(get_free_params(kernel))` long, in the
#'   same order. Each entry is the same shape as `evaluate(kernel, x1, x2)`.
#' @export
#' @examples
#' k <- se_kernel(length_scale = 2)
#' x <- c(0, 1, 3)
#' kernel_grad_free(k, x)[[1]]
#'
#' # Insensitive to a bare-name collision that would make kupdate() ambiguous:
#' k2 <- se_kernel(length_scale = 1) + se_kernel(length_scale = 3)
#' length(kernel_grad_free(k2, x)) == length(get_free_params(k2))
kernel_grad_free <- function(kernel, x1, x2 = NULL) {
  evaluate(kernel, x1, x2)

  raw <- kernel_grad_raw(kernel, x1, x2)
  params <- Filter(function(p) !isTRUE(p$frozen), flatten_params(kernel))

  out <- vector("list", 0)
  # raw and params stay aligned entry-for-entry for the same reason
  # kernel_grad()'s do (see its own comment): both walk own-fields-then-
  # subkernels in the same order.
  for (i in seq_along(params)) {
    entry <- raw[[i]]
    parametrisation <- params[[i]]$parametrisation
    wrapped <- params[[i]]$wrapped
    if (is.list(entry)) {
      for (k in seq_along(entry)) {
        out[[length(out) + 1]] <- entry[[k]] * unwrap_grad(parametrisation, wrapped[k])
      }
    } else {
      out[[length(out) + 1]] <- entry * unwrap_grad(parametrisation, wrapped)
    }
  }
  out
}
