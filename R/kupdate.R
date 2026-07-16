#' Immutably update a kernel's hyperparameters
#'
#' @description
#' Returns a copy of `kernel` with the named hyperparameter(s) updated to
#' new values, applied wherever they exist in the (possibly nested) kernel
#' tree. Validation of the new value happens inside `wrap()` of the relevant
#' [parametrisation] -- not duplicated here.
#'
#' Named `kupdate()` rather than `replace()` because `base::replace()`
#' already exists in R and is not an S3 generic (it has no `UseMethod()`
#' call), so overriding it would mean shadowing a base function -- avoided
#' deliberately.
#'
#' `kupdate()` is a single generic that works for every kernel, base or
#' composite, by walking `kernel$subkernels` (see [flatten_params()]),
#' rather than requiring each kernel type to implement its own update
#' logic.
#'
#' If a named argument does not match *any* hyperparameter, anywhere in the
#' kernel tree, `kupdate()` raises an error rather than silently ignoring it
#' (catches typos early). Likewise, a hyperparameter [frozen][freeze()] in
#' *every* occurrence where its name appears cannot be updated -- `kupdate()`
#' errors rather than silently no-op'ing (call [unfreeze()] first if that is
#' genuinely intended). If a name has a mix of frozen and non-frozen
#' occurrences in the tree (e.g. [ard_kernel()]'s internal, deliberately
#' frozen `length_scale` alongside an unrelated, trainable `length_scale`
#' elsewhere in a larger composite kernel), only the non-frozen occurrences
#' are updated -- frozen ones are left untouched, never silently overwritten.
#'
#' @param kernel A kernel object.
#' @param ... Named hyperparameter updates, e.g. `length_scale = 2`.
#'
#' @return A new kernel object with the requested updates applied.
#' @export
#' @examples
#' k <- se_kernel(length_scale = 1)
#' k2 <- kupdate(k, length_scale = 5)
#' get_params(k2)
kupdate <- function(kernel, ...) {
  updates <- list(...)
  if (length(updates) == 0) {
    return(kernel)
  }

  params <- flatten_params(kernel)
  available <- names(params)
  missing <- setdiff(names(updates), available)
  if (length(missing) > 0) {
    stop(
      "No such hyperparameter(s) in this kernel: ", paste(missing, collapse = ", "),
      ". Available: ", paste(unique(available), collapse = ", "), ".",
      call. = FALSE
    )
  }

  trainable <- names(Filter(function(p) !isTRUE(p$frozen), params))
  all_frozen <- setdiff(names(updates), trainable)
  if (length(all_frozen) > 0) {
    stop(
      "Hyperparameter(s) frozen in every occurrence, cannot kupdate(): ",
      paste(all_frozen, collapse = ", "),
      ". Call unfreeze() first if this is intentional.",
      call. = FALSE
    )
  }

  .kupdate_recurse(kernel, updates)
}

.kupdate_recurse <- function(kernel, updates) {
  for (name in names(kernel)) {
    field <- kernel[[name]]
    if (inherits(field, "wrapped_param") && name %in% names(updates) && !isTRUE(field$frozen)) {
      new_value <- updates[[name]]
      current_len <- length(field$wrapped)
      # A length mismatch that stays > 1 (e.g. 3 -> 2) is still caught, just
      # later, when evaluate.batch_kernel()/block_kernel()/... re-validates
      # its axis contract (.resolve_axes()) -- a clear error either way. But
      # collapsing to exactly 1 is different: .resolve_axes() treats a
      # length-1 value as a legitimate *shared* hyperparameter (that's how
      # one is declared in the first place), so it can never distinguish
      # "always meant to be shared" from "used to vary per batch/block
      # element, just silently stopped" -- once collapsed, this kupdate()
      # call is the only place that still knows which case it is. Caught
      # here, or not at all.
      if (current_len > 1 && length(new_value) == 1) {
        stop(
          "`", name, "` currently has ", current_len, " value(s) (e.g. a batch/block ",
          "axis) but kupdate() was given a single value -- this would silently turn it ",
          "into a shared value instead, discarding the per-element structure with no ",
          "error. Give exactly ", current_len, " replacement value(s) (use rep(x, ",
          current_len, ") if you intend the same value everywhere).",
          call. = FALSE
        )
      }
      kernel[[name]]$wrapped <- wrap(field$parametrisation, new_value)
    }
  }
  if (!is.null(kernel$subkernels)) {
    kernel$subkernels <- lapply(kernel$subkernels, .kupdate_recurse, updates = updates)
  }
  kernel
}

#' Set a kernel's free hyperparameters from the unconstrained (wrapped) space
#'
#' @description
#' The write-side counterpart of [get_free_params()]: assigns `values`
#' positionally, in the same tree-walk order [get_free_params()] reads in,
#' directly to each non-frozen hyperparameter's `$wrapped` field -- no
#' [wrap()] is applied, so `set_free_params(k, get_free_params(k))`
#' round-trips bit-exactly back to `k`. `values` must have exactly
#' [get_free_params()]'s length; a multi-valued hyperparameter (e.g.
#' [ard_kernel()]'s `length_scales`) consumes as many consecutive elements
#' as it has values.
#'
#' @section Numerical safety:
#' Only finiteness is checked here -- there is no upper bound on `values`,
#' by design (the unconstrained space is, in general, unbounded). This
#' matters for an optimiser driving `values` directly (e.g. `optim()` over
#' [get_free_params()]'s output): a single overly large step can produce a
#' wrapped value whose unwrapped (constrained) counterpart silently
#' overflows to an astronomically large but still-finite number -- e.g. a
#' [log_exp_parametrisation()] value of `700` unwraps to `exp(700) ~ 1e304`,
#' close to `double`'s range limit but not `Inf` -- which [evaluate()] and
#' [robust_chol()] will happily compute with, no warning, no error. Only a
#' large enough step to actually reach `Inf`/`NaN` is caught (downstream, by
#' [chol()]/[evaluate()] themselves, not here). If driving `values` with an
#' unbounded optimiser, consider a method that accepts box constraints
#' (e.g. `optim(..., method = "L-BFGS-B", lower = , upper = )`) to keep
#' steps within a sane range for your problem.
#'
#' @param kernel A kernel object.
#' @param values A finite numeric vector, `length(get_free_params(kernel))`
#'   long.
#' @return A new kernel object with the free hyperparameters updated.
#' @export
#' @examples
#' k <- se_kernel(length_scale = 1)
#' free <- get_free_params(k)
#' k2 <- set_free_params(k, free + 1) # a step in unconstrained space
#' get_trainable_params(k2)
set_free_params <- function(kernel, values) {
  if (!is.numeric(values) || anyNA(values) || any(!is.finite(values))) {
    stop("`values` must be a finite numeric vector.", call. = FALSE)
  }
  n_free <- length(get_free_params(kernel))
  if (length(values) != n_free) {
    stop(
      "`values` must have length ", n_free, " (the number of free/non-frozen ",
      "hyperparameter values in this kernel), got ", length(values), ".",
      call. = FALSE
    )
  }
  .set_free_params_recurse(kernel, values, 1L)$kernel
}

#' @keywords internal
.set_free_params_recurse <- function(kernel, values, pos) {
  for (name in names(kernel)) {
    field <- kernel[[name]]
    if (inherits(field, "wrapped_param") && !isTRUE(field$frozen)) {
      n <- length(field$wrapped)
      kernel[[name]]$wrapped <- values[pos:(pos + n - 1)]
      pos <- pos + n
    }
  }
  if (!is.null(kernel$subkernels)) {
    for (i in seq_along(kernel$subkernels)) {
      sub_result <- .set_free_params_recurse(kernel$subkernels[[i]], values, pos)
      kernel$subkernels[[i]] <- sub_result$kernel
      pos <- sub_result$pos
    }
  }
  list(kernel = kernel, pos = pos)
}
