#' Freeze or unfreeze hyperparameters
#'
#' @description
#' Marks named hyperparameters as excluded (`freeze()`) or included
#' (`unfreeze()`) from [get_trainable_params()], wherever they exist in the
#' (possibly nested) kernel tree. The `frozen` flag lives directly on the
#' relevant [wrapped_param()] object.
#'
#' `frozen` requires no parallel structure to keep in sync: it is metadata
#' carried directly by the parameter itself.
#'
#' @param kernel A kernel object.
#' @param ... Names of hyperparameters to freeze/unfreeze, e.g.
#'   `"length_scale"`.
#'
#' @return A new kernel object with the requested `frozen` flags applied.
#' @export
#' @examples
#' k <- se_kernel(length_scale = 1)
#' k_frozen <- freeze(k, "length_scale")
#' get_trainable_params(k_frozen)
freeze <- function(kernel, ...) {
  .set_frozen(kernel, c(...), TRUE)
}

#' @rdname freeze
#' @export
unfreeze <- function(kernel, ...) {
  .set_frozen(kernel, c(...), FALSE)
}

.set_frozen <- function(kernel, names_to_set, frozen) {
  if (length(names_to_set) == 0) {
    return(kernel)
  }
  available <- names(flatten_params(kernel))
  missing <- setdiff(names_to_set, available)
  if (length(missing) > 0) {
    stop(
      "No such hyperparameter(s) in this kernel: ", paste(missing, collapse = ", "),
      ". Available: ", paste(unique(available), collapse = ", "), ".",
      call. = FALSE
    )
  }
  .set_frozen_recurse(kernel, names_to_set, frozen)
}

.set_frozen_recurse <- function(kernel, names_to_set, frozen) {
  for (name in names(kernel)) {
    field <- kernel[[name]]
    if (inherits(field, "wrapped_param") && name %in% names_to_set) {
      kernel[[name]]$frozen <- frozen
    }
  }
  if (!is.null(kernel$subkernels)) {
    kernel$subkernels <- lapply(
      kernel$subkernels, .set_frozen_recurse,
      names_to_set = names_to_set, frozen = frozen
    )
  }
  kernel
}
