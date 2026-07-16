# Operator kernels: binary compositions, stored as $subkernels (length 2),
# with no engine of their own -- they override evaluate() directly and
# recombine the evaluate() results of their sub-kernels (layer 4 of the
# architecture, see dev/ameliorations/design/01-architecture.Rmd). This is also why no new
# code is needed in flatten_params()/kupdate()/freeze() for them: those
# already recurse through $subkernels generically.
#
# Per the user's explicit decision, `-` (subtraction)/unary negation are not
# implemented -- no neg_kernel(). Only `+` and `*` are supported operators.

#' Coerce a value to a kernel
#'
#' @description
#' Used by [sum_kernel()]/[product_kernel()] (and the `+`/`*` operators), and
#' by every wrapper taking an `inner` kernel, to let a bare number stand in
#' for a constant kernel. Only a `kernel` object or a single finite numeric
#' value is accepted -- anything else (`NULL`, `NA`, a vector, a character
#' string...) is rejected here, at construction, rather than left to fail
#' later with a confusing low-level error once the malformed object reaches
#' [evaluate()].
#'
#' @param x A kernel object, or a single finite numeric value.
#' @return A kernel object.
#' @keywords internal
as_kernel <- function(x) {
  if (inherits(x, "kernel")) {
    return(x)
  }
  if (!is.numeric(x) || length(x) != 1 || !is.finite(x)) {
    stop(
      "Expected a `kernel` object or a single finite numeric value ",
      "(coerced to a constant_kernel()), got ",
      if (is.null(x)) {
        "NULL"
      } else {
        paste0("`", class(x)[1], "`", if (length(x) > 1) paste0(" of length ", length(x)) else "")
      },
      ".",
      call. = FALSE
    )
  }
  constant_kernel(value = x)
}

#' Sum of two kernels
#'
#' @description
#' \deqn{k(x_1, x_2) = k_1(x_1, x_2) + k_2(x_1, x_2)}
#'
#' A bare number is automatically coerced to a [constant_kernel()] (see
#' [as_kernel()]). Usually constructed via the `+` operator (`k1 + k2`)
#' rather than called directly.
#'
#' @param left,right Kernel objects, or numeric values (coerced via
#'   [as_kernel()]).
#'
#' @return An object of class `sum_kernel`/`operator_kernel`/`kernel`.
#' @export
#' @examples
#' k <- se_kernel(length_scale = 1) + white_noise_kernel(noise = 0.1)
#' evaluate(k, c(0, 1, 2))
sum_kernel <- function(left, right) {
  kernel <- structure(
    list(subkernels = list(as_kernel(left), as_kernel(right))),
    class = c("sum_kernel", "operator_kernel", "kernel")
  )
  validate_kernel(kernel)
  kernel
}

#' @export
evaluate.sum_kernel <- function(kernel, x1, x2 = NULL, ...) {
  left <- evaluate(kernel$subkernels[[1]], x1, x2, ...)
  right <- evaluate(kernel$subkernels[[2]], x1, x2, ...)
  .check_combinable_shapes(left, right, "sum_kernel")
  left + right
}

#' @export
format.sum_kernel <- function(x, ...) {
  paste(format(x$subkernels[[1]]), "+", format(x$subkernels[[2]]))
}

#' Product of two kernels
#'
#' @description
#' \deqn{k(x_1, x_2) = k_1(x_1, x_2) \cdot k_2(x_1, x_2)}
#'
#' A bare number is automatically coerced to a [constant_kernel()] (see
#' [as_kernel()]). Usually constructed via the `*` operator (`k1 * k2`)
#' rather than called directly. A common use of `product_kernel()` is to
#' scale a kernel with a [variance_kernel()].
#'
#' @param left,right Kernel objects, or numeric values (coerced via
#'   [as_kernel()]).
#'
#' @return An object of class `product_kernel`/`operator_kernel`/`kernel`.
#' @export
#' @examples
#' k <- variance_kernel(variance = 2) * se_kernel(length_scale = 1)
#' evaluate(k, c(0, 1, 2))
product_kernel <- function(left, right) {
  kernel <- structure(
    list(subkernels = list(as_kernel(left), as_kernel(right))),
    class = c("product_kernel", "operator_kernel", "kernel")
  )
  validate_kernel(kernel)
  kernel
}

#' @export
evaluate.product_kernel <- function(kernel, x1, x2 = NULL, ...) {
  left <- evaluate(kernel$subkernels[[1]], x1, x2, ...)
  right <- evaluate(kernel$subkernels[[2]], x1, x2, ...)
  .check_combinable_shapes(left, right, "product_kernel")
  left * right
}

#' Check that two sub-kernel results can be combined elementwise
#'
#' @description
#' `sum_kernel`/`product_kernel` combine their two sub-kernels' `evaluate()`
#' results with `+`/`*`. If one sub-kernel uses [dense_engine()] (a matrix
#' result) and the other [diagonal_engine()] (a vector result), R's ordinary
#' vector-recycling rules would combine them anyway -- silently producing a
#' result that adds the vector to every row of the matrix, which is not a
#' meaningful covariance for any GP interpretation (and, unlike a genuine
#' shape mismatch, does not always error: with a square, symmetric operand
#' it can even pass [check_kernel_contract()]'s symmetry/PSD checks by
#' accident). Caught here, once, rather than left to silently corrupt a
#' composite kernel's covariance.
#'
#' @param left,right The two sub-kernels' `evaluate()` results.
#' @param kernel_class The composite kernel's class, for the error message.
#' @keywords internal
.check_combinable_shapes <- function(left, right, kernel_class) {
  left_is_matrix <- is.matrix(left)
  right_is_matrix <- is.matrix(right)
  if (left_is_matrix != right_is_matrix) {
    stop(
      "`", kernel_class, "` cannot combine a sub-kernel using dense_engine() ",
      "(a matrix result) with one using diagonal_engine() (a vector result) -- ",
      "give both sub-kernels the same engine.",
      call. = FALSE
    )
  }
  if (left_is_matrix && !identical(dim(left), dim(right))) {
    stop(
      "`", kernel_class, "` sub-kernels produced incompatible shapes: ",
      paste(dim(left), collapse = "x"), " vs ", paste(dim(right), collapse = "x"), ".",
      call. = FALSE
    )
  }
  if (!left_is_matrix && length(left) != length(right)) {
    stop(
      "`", kernel_class, "` sub-kernels produced vectors of different length: ",
      length(left), " vs ", length(right), ".",
      call. = FALSE
    )
  }
}

#' @export
format.product_kernel <- function(x, ...) {
  # Any sum_kernel operand is always parenthesised when formatting a
  # product, so `(k1 + k2) * k3` prints unambiguously rather than risking
  # a string that could be misread as k1 + (k2 * k3).
  fmt_side <- function(k) {
    s <- format(k)
    if (inherits(k, "sum_kernel")) paste0("(", s, ")") else s
  }
  paste(fmt_side(x$subkernels[[1]]), "*", fmt_side(x$subkernels[[2]]))
}

#' Arithmetic operators for kernels
#'
#' @description
#' `k1 + k2` builds a [sum_kernel()]; `k1 * k2` builds a [product_kernel()].
#' Either side may be a bare number, automatically coerced to a
#' [constant_kernel()] (see [as_kernel()]).
#'
#' `-` (subtraction) and unary negation are deliberately not implemented --
#' there is no `neg_kernel()`. Subtracting one covariance kernel from
#' another rarely yields a valid (positive semi-definite) covariance
#' function, so it is omitted rather than silently producing an invalid
#' kernel.
#'
#' @param e1,e2 Kernel objects, or numeric values.
#' @return A new kernel object.
#' @name kernel_arithmetic
#' @export
#' @examples
#' k <- se_kernel(length_scale = 1) + 0.5
#' evaluate(k, c(0, 1, 2))
Ops.kernel <- function(e1, e2) {
  if (.Generic == "+") {
    return(sum_kernel(e1, e2))
  }
  if (.Generic == "*") {
    return(product_kernel(e1, e2))
  }
  stop(
    "Operator `", .Generic, "` is not defined for kernel objects ",
    "(only `+` and `*` are supported -- see ?kernel_arithmetic).",
    call. = FALSE
  )
}
