#' @export
format.kernel <- function(x, ...) {
  params <- Filter(function(field) inherits(field, "wrapped_param"), x)
  parts <- vapply(
    names(params),
    function(name) paste0(name, " = ", format(params[[name]])),
    character(1)
  )
  paste0(class(x)[[1]], "(", paste(parts, collapse = ", "), ")")
}

#' @export
print.kernel <- function(x, ...) {
  cat(format(x), "\n", sep = "")
  invisible(x)
}

#' Print a kernel as an indented tree of classes and hyperparameters
#'
#' @description
#' `format()`/`print()` give a single-line, composable string (e.g.
#' `"k1 + k2"`) meant to stay embeddable in error messages and plot titles.
#' `kernel_tree()` is a complementary, more verbose view: one line per node
#' of the `$subkernels` tree, indented by nesting depth, each with its own
#' class name and hyperparameters -- values and [frozen][freeze()] status
#' (via the same `format.wrapped_param()` used everywhere else, see
#' [wrapped_param()]). Useful to inspect a deeply nested composite kernel at
#' a glance.
#'
#' @param kernel A kernel object.
#' @param indent Internal use (recursion depth).
#' @return `kernel`, invisibly (called for its printed output).
#' @export
#' @examples
#' k <- (variance_kernel(0.4) * se_kernel(length_scale = 2)) +
#'   freeze(white_noise_kernel(noise = 0.05), "noise")
#' kernel_tree(k)
kernel_tree <- function(kernel, indent = 0) {
  params <- Filter(function(field) inherits(field, "wrapped_param"), kernel)
  hp_txt <- if (length(params) == 0) {
    ""
  } else {
    parts <- vapply(
      names(params),
      function(name) paste0(name, " = ", format(params[[name]])),
      character(1)
    )
    paste0(" (", paste(parts, collapse = ", "), ")")
  }
  cat(strrep("  ", indent), class(kernel)[1], hp_txt, "\n", sep = "")
  for (sub in kernel$subkernels) {
    kernel_tree(sub, indent + 1)
  }
  invisible(kernel)
}
