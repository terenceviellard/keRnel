# Dot-product kernels share the same "shape() of a scalar transform" pattern
# as stationary kernels (R/stationary-kernels.R), but the scalar in question
# is a dot product rather than a distance: pairwise(kernel, x1, x2) is always
# shape(kernel, distance_function(x1, x2)), with distance_function defaulting
# to dot_product(). This lets pairwise_matrix.dot_product_kernel()
# (R/distance-matrix.R) reuse the vectorised dot_product fast path
# (X1 %*% t(X2)) the same way pairwise_matrix.stationary_kernel() reuses the
# vectorised distance fast paths.

#' Linear kernel
#'
#' @description
#' \deqn{k(x_1, x_2) = \mathrm{slope\_var} \cdot \langle x_1, x_2 \rangle}
#'
#' Note: samples/posteriors from this kernel always cross the point
#' `(0, 0)`. To allow other crossing points, either use [affine_kernel()]
#' with a fixed `offset`, or add a [constant_kernel()] to a
#' `linear_kernel()`, where the constant value represents the variance at
#' the crossing point.
#'
#' @param slope_var A strictly positive numeric value. Controls the slope.
#' @param slope_var_parametrisation A [parametrisation] for `slope_var`.
#'   Defaults to [log_exp_parametrisation()].
#' @param distance_function A distance function, or a shortcut name (see
#'   [resolve_distance_function()]). Defaults to [dot_product()].
#' @param engine A computation engine. Defaults to [dense_engine()].
#'
#' @return An object of class `linear_kernel`/`dot_product_kernel`/`kernel`.
#' @export
#' @examples
#' k <- linear_kernel(slope_var = 1)
#' evaluate(k, c(0, 1, 2))
linear_kernel <- function(slope_var,
                           slope_var_parametrisation = log_exp_parametrisation(),
                           distance_function = dot_product,
                           engine = dense_engine()) {
  kernel <- structure(
    list(
      slope_var = wrapped_param(slope_var, slope_var_parametrisation),
      distance_function = resolve_distance_function(distance_function),
      engine = engine
    ),
    class = c("linear_kernel", "dot_product_kernel", "kernel")
  )
  validate_kernel(kernel)
  kernel
}

#' @export
pairwise.linear_kernel <- function(kernel, x1, x2) {
  shape(kernel, kernel$distance_function(x1, x2))
}

#' @export
shape.linear_kernel <- function(kernel, d) {
  unwrap_param(kernel$slope_var) * d
}

#' @export
shape_grad.linear_kernel <- function(kernel, d) {
  list(slope_var = d)
}

#' @export
shape_grad_d.linear_kernel <- function(kernel, d) {
  0 * d + unwrap_param(kernel$slope_var)
}

#' Affine kernel
#'
#' @description
#' \deqn{k(x_1, x_2) = \mathrm{slope\_var} \cdot \langle x_1 - \mathrm{offset}, x_2 - \mathrm{offset} \rangle}
#'
#' Note: samples/posteriors from this kernel always cross the point
#' `(offset, 0)`. To add uncertainty about the crossing point, add a
#' [constant_kernel()] to an `affine_kernel()`.
#'
#' @param slope_var A strictly positive numeric value. Controls the slope.
#' @param offset A numeric value (any real value): determines the crossing
#'   point. Not optimisable via [get_params()]: it is a plain field, not a
#'   hyperparameter.
#' @param slope_var_parametrisation A [parametrisation] for `slope_var`.
#'   Defaults to [log_exp_parametrisation()].
#' @param distance_function A distance function, or a shortcut name (see
#'   [resolve_distance_function()]). Defaults to [dot_product()].
#' @param engine A computation engine. Defaults to [dense_engine()].
#'
#' @return An object of class `affine_kernel`/`dot_product_kernel`/`kernel`.
#' @export
#' @examples
#' k <- affine_kernel(slope_var = 1, offset = 0.5)
#' evaluate(k, c(0, 1, 2))
affine_kernel <- function(slope_var,
                           offset,
                           slope_var_parametrisation = log_exp_parametrisation(),
                           distance_function = dot_product,
                           engine = dense_engine()) {
  kernel <- structure(
    list(
      slope_var = wrapped_param(slope_var, slope_var_parametrisation),
      offset = offset,
      distance_function = resolve_distance_function(distance_function),
      engine = engine
    ),
    class = c("affine_kernel", "dot_product_kernel", "kernel")
  )
  validate_kernel(kernel)
  kernel
}

#' @export
pairwise.affine_kernel <- function(kernel, x1, x2) {
  shape(kernel, kernel$distance_function(x1 - kernel$offset, x2 - kernel$offset))
}

#' @export
shape.affine_kernel <- function(kernel, d) {
  unwrap_param(kernel$slope_var) * d
}

#' @export
shape_grad.affine_kernel <- function(kernel, d) {
  list(slope_var = d)
}

#' @export
shape_grad_d.affine_kernel <- function(kernel, d) {
  0 * d + unwrap_param(kernel$slope_var)
}

#' @keywords internal
#' @exportS3Method
pairwise_matrix.affine_kernel <- function(kernel, X1, X2) {
  X1c <- sweep(X1, 2, kernel$offset, "-")
  X2c <- sweep(X2, 2, kernel$offset, "-")
  shape(kernel, distance_matrix(kernel$distance_function, X1c, X2c))
}

#' @keywords internal
#' @exportS3Method
pairwise_diag.affine_kernel <- function(kernel, X1, X2) {
  X1c <- sweep(X1, 2, kernel$offset, "-")
  X2c <- sweep(X2, 2, kernel$offset, "-")
  shape(kernel, diag_distance(kernel$distance_function, X1c, X2c))
}

# affine_kernel() centers by `offset` before computing the distance (see
# pairwise_matrix.affine_kernel() above), so it needs its own gradient
# override too -- the pairwise_matrix_grad.dot_product_kernel() fallback
# would otherwise differentiate the *uncentered* distance instead.

#' @keywords internal
#' @exportS3Method
pairwise_matrix_grad.affine_kernel <- function(kernel, X1, X2) {
  X1c <- sweep(X1, 2, kernel$offset, "-")
  X2c <- sweep(X2, 2, kernel$offset, "-")
  shape_grad(kernel, distance_matrix(kernel$distance_function, X1c, X2c))
}

#' @keywords internal
#' @exportS3Method
pairwise_diag_grad.affine_kernel <- function(kernel, X1, X2) {
  X1c <- sweep(X1, 2, kernel$offset, "-")
  X2c <- sweep(X2, 2, kernel$offset, "-")
  shape_grad(kernel, diag_distance(kernel$distance_function, X1c, X2c))
}

#' Polynomial kernel
#'
#' @description
#' \deqn{k(x_1, x_2) = (\gamma \cdot \langle x_1, x_2 \rangle + c)^{\mathrm{degree}}}
#'
#' @param degree A positive integer. Not optimisable via [get_params()] or
#'   [kupdate()] -- create a new kernel instance to change it.
#' @param gamma A strictly positive numeric value.
#' @param constant A numeric value (any real value): the independent term.
#'   Not optimisable via [get_params()].
#' @param gamma_parametrisation A [parametrisation] for `gamma`. Defaults to
#'   [log_exp_parametrisation()].
#' @param distance_function A distance function, or a shortcut name (see
#'   [resolve_distance_function()]). Defaults to [dot_product()].
#' @param engine A computation engine. Defaults to [dense_engine()].
#'
#' @return An object of class `polynomial_kernel`/`dot_product_kernel`/
#'   `kernel`.
#' @export
#' @examples
#' k <- polynomial_kernel(degree = 2, gamma = 1, constant = 1)
#' evaluate(k, c(0, 1, 2))
polynomial_kernel <- function(degree,
                               gamma,
                               constant,
                               gamma_parametrisation = log_exp_parametrisation(),
                               distance_function = dot_product,
                               engine = dense_engine()) {
  if (!is.numeric(degree) || length(degree) != 1 || degree <= 0 || degree != round(degree)) {
    stop("`degree` must be a single positive integer.", call. = FALSE)
  }

  kernel <- structure(
    list(
      degree = as.integer(degree),
      gamma = wrapped_param(gamma, gamma_parametrisation),
      constant = constant,
      distance_function = resolve_distance_function(distance_function),
      engine = engine
    ),
    class = c("polynomial_kernel", "dot_product_kernel", "kernel")
  )
  validate_kernel(kernel)
  kernel
}

#' @export
pairwise.polynomial_kernel <- function(kernel, x1, x2) {
  shape(kernel, kernel$distance_function(x1, x2))
}

#' @export
shape.polynomial_kernel <- function(kernel, d) {
  gamma <- unwrap_param(kernel$gamma)
  (gamma * d + kernel$constant)^kernel$degree
}

#' @export
shape_grad.polynomial_kernel <- function(kernel, d) {
  gamma <- unwrap_param(kernel$gamma)
  base <- gamma * d + kernel$constant
  list(gamma = kernel$degree * d * base^(kernel$degree - 1))
}

#' @export
shape_grad_d.polynomial_kernel <- function(kernel, d) {
  gamma <- unwrap_param(kernel$gamma)
  base <- gamma * d + kernel$constant
  kernel$degree * gamma * base^(kernel$degree - 1)
}

#' Sigmoid kernel
#'
#' @description
#' \deqn{k(x_1, x_2) = \tanh(\alpha \cdot \langle x_1, x_2 \rangle + c)}
#'
#' @param alpha A strictly positive numeric value.
#' @param constant A numeric value (any real value).
#' @param alpha_parametrisation A [parametrisation] for `alpha`. Defaults to
#'   [log_exp_parametrisation()].
#' @param distance_function A distance function, or a shortcut name (see
#'   [resolve_distance_function()]). Defaults to [dot_product()].
#' @param engine A computation engine. Defaults to [dense_engine()].
#'
#' @return An object of class `sigmoid_kernel`/`dot_product_kernel`/`kernel`.
#' @export
#' @examples
#' k <- sigmoid_kernel(alpha = 1, constant = 0)
#' evaluate(k, c(0, 1, 2))
sigmoid_kernel <- function(alpha,
                            constant,
                            alpha_parametrisation = log_exp_parametrisation(),
                            distance_function = dot_product,
                            engine = dense_engine()) {
  kernel <- structure(
    list(
      alpha = wrapped_param(alpha, alpha_parametrisation),
      constant = constant,
      distance_function = resolve_distance_function(distance_function),
      engine = engine
    ),
    class = c("sigmoid_kernel", "dot_product_kernel", "kernel")
  )
  validate_kernel(kernel)
  kernel
}

#' @export
pairwise.sigmoid_kernel <- function(kernel, x1, x2) {
  shape(kernel, kernel$distance_function(x1, x2))
}

#' @export
shape.sigmoid_kernel <- function(kernel, d) {
  tanh(unwrap_param(kernel$alpha) * d + kernel$constant)
}

#' @export
shape_grad.sigmoid_kernel <- function(kernel, d) {
  s <- shape(kernel, d)
  list(alpha = d * (1 - s^2))
}

#' @export
shape_grad_d.sigmoid_kernel <- function(kernel, d) {
  s <- shape(kernel, d)
  unwrap_param(kernel$alpha) * (1 - s^2)
}
