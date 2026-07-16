#' Constant kernel
#'
#' @description
#' Returns the same constant value regardless of the inputs:
#'
#' \deqn{k(x_1, x_2) = \mathrm{value}}
#'
#' Unlike [variance_kernel()], `value` is unconstrained: it may be any
#' finite real number (negative included). It is still stored as a
#' [wrapped_param()] (with [identity_parametrisation()] by default), for
#' architectural consistency with every other hyperparameter -- so `value`
#' is visible to [get_params()]/[kupdate()]/[freeze()] like any other
#' optimisable hyperparameter, and validated (finite, numeric) at
#' construction and at every [kupdate()].
#'
#' @param value A finite numeric value. Defaults to `1`.
#' @param value_parametrisation A [parametrisation] for `value`. Defaults to
#'   [identity_parametrisation()].
#' @param engine A computation engine. Defaults to [dense_engine()].
#'
#' @return An object of class `constant_kernel`/`kernel`.
#' @export
#' @examples
#' k <- constant_kernel(value = 2)
#' evaluate(k, c(0, 1, 2))
constant_kernel <- function(value = 1,
                             value_parametrisation = identity_parametrisation(),
                             engine = dense_engine()) {
  kernel <- structure(
    list(value = wrapped_param(value, value_parametrisation), engine = engine),
    class = c("constant_kernel", "kernel")
  )
  validate_kernel(kernel)
  kernel
}

#' @export
pairwise.constant_kernel <- function(kernel, x1, x2) {
  unwrap_param(kernel$value)
}

#' @keywords internal
#' @exportS3Method
pairwise_matrix.constant_kernel <- function(kernel, X1, X2) {
  matrix(unwrap_param(kernel$value), nrow = nrow(X1), ncol = nrow(X2))
}

#' Variance kernel
#'
#' @description
#' Returns the same strictly positive constant value everywhere. Typically
#' multiplied with another kernel (via `*`, see [kernel_arithmetic]) to
#' scale its output.
#'
#' \deqn{k(x_1, x_2) = \mathrm{variance}}
#'
#' @param variance A strictly positive numeric value. Defaults to `1`.
#' @param variance_parametrisation A [parametrisation] for `variance`.
#'   Defaults to [log_exp_parametrisation()].
#' @param engine A computation engine. Defaults to [dense_engine()].
#'
#' @return An object of class `variance_kernel`/`kernel`.
#' @export
#' @examples
#' k <- variance_kernel(variance = 2)
#' evaluate(k, c(0, 1, 2))
variance_kernel <- function(variance = 1,
                             variance_parametrisation = log_exp_parametrisation(),
                             engine = dense_engine()) {
  kernel <- structure(
    list(
      variance = wrapped_param(variance, variance_parametrisation),
      engine = engine
    ),
    class = c("variance_kernel", "kernel")
  )
  validate_kernel(kernel)
  kernel
}

#' @export
pairwise.variance_kernel <- function(kernel, x1, x2) {
  unwrap_param(kernel$variance)
}

#' @keywords internal
#' @exportS3Method
pairwise_matrix.variance_kernel <- function(kernel, X1, X2) {
  matrix(unwrap_param(kernel$variance), nrow = nrow(X1), ncol = nrow(X2))
}
