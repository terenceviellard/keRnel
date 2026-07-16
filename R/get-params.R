# The single recursive tree-walk that the whole architecture relies on (see
# dev/ameliorations/design/01-architecture.Rmd, section on $subkernels). flatten_params() is
# its canonical, read-only form, reused (with different logic at the leaves)
# by get_params()/get_trainable_params() here, by kupdate() (R/kupdate.R),
# by freeze()/unfreeze() (R/freeze.R), and -- once batching exists --
# by slice_kernel().

#' Recursively collect every `wrapped_param` leaf of a kernel
#'
#' @description
#' Walks `kernel`'s own fields and, if present, recurses into
#' `kernel$subkernels`. Returns a flat named list of every `wrapped_param`
#' object found, named after their field name at each level (sub-kernel
#' field names are not currently disambiguated by depth/position -- with a
#' single level of nesting this is unambiguous; revisit once kernels with
#' deeper nesting and same-named parameters are implemented).
#'
#' @param kernel A kernel object.
#' @return A named list of `wrapped_param` objects.
#' @keywords internal
flatten_params <- function(kernel) {
  own <- Filter(function(field) inherits(field, "wrapped_param"), kernel)
  nested <- if (!is.null(kernel$subkernels)) {
    unlist(lapply(kernel$subkernels, flatten_params), recursive = FALSE)
  } else {
    list()
  }
  c(own, nested)
}

#' Like `flatten_params()`, but also records each `wrapped_param`'s owning
#' kernel's class alongside it.
#'
#' @description
#' Used only to build human-readable, class-qualified names for
#' [get_params()]/[get_trainable_params()] when a hyperparameter name is
#' shared by sibling sub-kernels (see [check_param_names()]) -- deliberately
#' **not** merged into [flatten_params()] itself, which is the shared
#' tree-walk skeleton reused by [kupdate()]/[freeze()]/[sample_hps()]/
#' `.resolve_axes()`; none of those need, or should be affected by, this
#' extra bookkeeping.
#'
#' @param kernel A kernel object.
#' @return A list with `params` (same as [flatten_params()]) and `classes`
#'   (a character vector, `class(owner)[1]` for each entry, same length and
#'   order as `params`).
#' @keywords internal
.flatten_params_with_class <- function(kernel) {
  own <- Filter(function(field) inherits(field, "wrapped_param"), kernel)
  result <- list(params = own, classes = rep(class(kernel)[1], length(own)))
  for (sub in kernel$subkernels) {
    nested <- .flatten_params_with_class(sub)
    result$params <- c(result$params, nested$params)
    result$classes <- c(result$classes, nested$classes)
  }
  result
}

#' Disambiguate hyperparameter names shared across sibling sub-kernels
#'
#' @description
#' A bare name (`length_scale`) is left untouched unless another occurrence
#' of the same name exists elsewhere in the tree. A colliding name is
#' qualified with its owning kernel's class (`length_scale (se_kernel)`,
#' `length_scale (periodic_kernel)`) -- the same class name already visible
#' in `format()`'s tree-shaped output, so the two stay consistent. If the
#' *same* class also repeats (e.g. `se_kernel(1) + se_kernel(2)`, where the
#' class alone cannot disambiguate), occurrences within that class are
#' additionally numbered (`length_scale (se_kernel #1)`,
#' `length_scale (se_kernel #2)`), scoped to the colliding group only --
#' never a global positional index over the whole tree.
#'
#' @param names Character vector of hyperparameter names.
#' @param classes Character vector, same length, the owning kernel class of
#'   each name.
#' @return `names`, with colliding entries qualified.
#' @keywords internal
.disambiguate_param_names <- function(names, classes) {
  out <- names
  for (nm in unique(names[duplicated(names)])) {
    idx <- which(names == nm)
    cls <- classes[idx]
    qualified <- paste0(nm, " (", cls, ")")
    repeats_class <- stats::ave(seq_along(idx), cls, FUN = length) > 1
    if (any(repeats_class)) {
      counter <- stats::ave(seq_along(idx), cls, FUN = seq_along)
      qualified[repeats_class] <- paste0(nm, " (", cls[repeats_class], " #", counter[repeats_class], ")")
    }
    out[idx] <- qualified
  }
  out
}

#' Get every hyperparameter of a kernel, in its natural (constrained) space
#'
#' @description
#' Flattens the kernel's hyperparameter tree and unwraps every value. The
#' result is suitable to inspect a kernel's current configuration. For the
#' values an optimiser should actually vary, see [get_trainable_params()],
#' which excludes parameters frozen via [freeze()].
#'
#' If a hyperparameter name is shared by more than one sub-kernel (e.g.
#' `se_kernel(1) + periodic_kernel(length_scale = 2, period = 3)`, both
#' having a `length_scale`), the shared entries are returned under
#' disambiguated names (`"length_scale (se_kernel)"`,
#' `"length_scale (periodic_kernel)"`) rather than as ambiguous duplicate
#' names in the returned vector -- see [.disambiguate_param_names()]. This
#' does not affect [kupdate()]/[freeze()], which continue to match by the
#' bare name, structurally, everywhere it occurs in the tree (their
#' documented, intentional "propagate to every occurrence" behaviour -- see
#' [kupdate()]). Call [check_param_names()] if you want to know, ahead of
#' time, whether a name collision in a kernel you built is likely
#' intentional (shared HP) or a typo.
#'
#' @param kernel A kernel object.
#' @return A named numeric vector.
#' @export
#' @examples
#' get_params(se_kernel(length_scale = 2))
get_params <- function(kernel) {
  fp <- .flatten_params_with_class(kernel)
  names(fp$params) <- .disambiguate_param_names(names(fp$params), fp$classes)
  unlist(lapply(fp$params, unwrap_param))
}

#' Get the non-frozen hyperparameters of a kernel
#'
#' @description
#' Same as [get_params()], but excludes any `wrapped_param` marked
#' `frozen = TRUE` (via [freeze()]). This is the vector an optimiser should
#' actually act on. Names are disambiguated the same way as [get_params()]
#' when shared across sub-kernels.
#'
#' @param kernel A kernel object.
#' @return A named numeric vector.
#' @export
#' @examples
#' k <- freeze(se_kernel(length_scale = 2), "length_scale")
#' get_params(k)            # length_scale still listed
#' get_trainable_params(k)  # length_scale excluded
get_trainable_params <- function(kernel) {
  fp <- .flatten_params_with_class(kernel)
  names(fp$params) <- .disambiguate_param_names(names(fp$params), fp$classes)
  trainable <- Filter(function(p) !isTRUE(p$frozen), fp$params)
  unlist(lapply(trainable, unwrap_param))
}

#' Get a kernel's free hyperparameters in the unconstrained (wrapped) space
#'
#' @description
#' Like [get_trainable_params()], but returns each non-frozen
#' hyperparameter's *wrapped* value (`$wrapped`, the unconstrained space a
#' [parametrisation] maps to/from) instead of its constrained value -- no
#' [unwrap_param()] is applied. Paired with [set_free_params()] for an
#' optimiser that works natively in unconstrained space (e.g. plain
#' gradient descent or `optim(method = "BFGS")` without box constraints):
#' round-tripping through `unwrap()`/`wrap()` on every iteration (the
#' constrained-space alternative) is not guaranteed bit-exact, since both
#' are floating-point transcendental functions (e.g. `log(exp(w))` need not
#' equal `w`) -- reading and writing `$wrapped` directly avoids that drift
#' entirely.
#'
#' The result is deliberately **unnamed and positional**, unlike
#' [get_trainable_params()]: two hyperparameters sharing a bare name at
#' different points in the tree (see [.disambiguate_param_names()]) are
#' still perfectly distinguishable here, by position, whereas a name-keyed
#' API would inherit the same ambiguity [kupdate()] has (it updates *every*
#' occurrence of a bare name at once). [set_free_params()] relies on this
#' same tree-walk order to write values back -- see its "Numerical safety"
#' section for a caveat about driving these values with an unbounded
#' optimiser.
#'
#' @param kernel A kernel object.
#' @return A plain (unnamed) numeric vector, in tree-walk order.
#' @export
#' @examples
#' k <- se_kernel(length_scale = 2)
#' get_free_params(k)         # wrap(log_exp_parametrisation(), 2) == log(2)
#' get_trainable_params(k)    # 2 (the constrained value)
get_free_params <- function(kernel) {
  params <- flatten_params(kernel)
  trainable <- Filter(function(p) !isTRUE(p$frozen), params)
  unlist(lapply(trainable, function(p) p$wrapped), use.names = FALSE)
}
