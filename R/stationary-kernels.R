#' Squared Exponential (RBF) kernel
#'
#' @description
#' The Squared Exponential kernel (also known as "RBF" or "Gaussian"):
#'
#' \deqn{k(x_1, x_2) = \exp\left(-\frac{d(x_1, x_2)}{2 \cdot \mathrm{length\_scale}^2}\right)}
#'
#' where \eqn{d} is `distance_function` (squared Euclidean by default).
#' `length_scale` is constrained to stay strictly positive.
#'
#' @param length_scale A strictly positive numeric value (or vector, for a
#'   batched kernel -- see [batch_kernel()]).
#' @param length_scale_parametrisation A [parametrisation] for
#'   `length_scale`. Defaults to [log_exp_parametrisation()].
#' @param distance_function A distance function, or a shortcut name (see
#'   [resolve_distance_function()]). Defaults to
#'   [squared_euclidean_distance()].
#' @param engine A computation engine. Defaults to [dense_engine()].
#'
#' @return An object of class `se_kernel`/`stationary_kernel`/`kernel`.
#' @export
#' @examples
#' k <- se_kernel(length_scale = 1)
#' evaluate(k, c(0, 1, 2))
se_kernel <- function(length_scale,
                       length_scale_parametrisation = log_exp_parametrisation(),
                       distance_function = squared_euclidean_distance,
                       engine = dense_engine()) {
  kernel <- structure(
    list(
      length_scale = wrapped_param(length_scale, length_scale_parametrisation),
      distance_function = resolve_distance_function(distance_function),
      engine = engine
    ),
    class = c("se_kernel", "stationary_kernel", "kernel")
  )
  validate_kernel(kernel)
  kernel
}

#' @rdname se_kernel
#' @export
rbf_kernel <- se_kernel

#' @export
pairwise.se_kernel <- function(kernel, x1, x2) {
  shape(kernel, kernel$distance_function(x1, x2))
}

#' @export
shape.se_kernel <- function(kernel, d) {
  length_scale <- unwrap_param(kernel$length_scale)
  exp(-0.5 * d / length_scale^2)
}

#' @export
shape_grad.se_kernel <- function(kernel, d) {
  length_scale <- unwrap_param(kernel$length_scale)
  s <- shape(kernel, d)
  list(length_scale = s * d / length_scale^3)
}

#' @export
shape_grad_d.se_kernel <- function(kernel, d) {
  length_scale <- unwrap_param(kernel$length_scale)
  -0.5 * shape(kernel, d) / length_scale^2
}

#' Matern kernels
#'
#' @description
#' The Matern family of stationary kernels, of smoothness 1/2 ("exponential"),
#' 3/2, and 5/2:
#'
#' \deqn{k_{1/2}(x_1, x_2) = \exp\left(-\frac{d}{\mathrm{length\_scale}}\right)}
#' \deqn{k_{3/2}(x_1, x_2) = \left(1 + \frac{\sqrt{3} d}{l}\right)\exp\left(-\frac{\sqrt{3} d}{l}\right)}
#' \deqn{k_{5/2}(x_1, x_2) = \left(1 + \frac{\sqrt{5} d}{l} + \frac{5 d^2}{3 l^2}\right)\exp\left(-\frac{\sqrt{5} d}{l}\right)}
#'
#' where \eqn{d} is `distance_function` (Euclidean by default) and
#' \eqn{l} is `length_scale`, constrained to stay strictly positive.
#'
#' @param length_scale A strictly positive numeric value.
#' @param length_scale_parametrisation A [parametrisation] for
#'   `length_scale`. Defaults to [log_exp_parametrisation()].
#' @param distance_function A distance function, or a shortcut name (see
#'   [resolve_distance_function()]). Defaults to [euclidean_distance()].
#' @param engine A computation engine. Defaults to [dense_engine()].
#'
#' @return An object of class `matern12_kernel`/`matern32_kernel`/
#'   `matern52_kernel`, `stationary_kernel`, `kernel`.
#' @name matern_kernels
#' @export
#' @examples
#' k <- matern32_kernel(length_scale = 1)
#' evaluate(k, c(0, 1, 2))
matern12_kernel <- function(length_scale,
                             length_scale_parametrisation = log_exp_parametrisation(),
                             distance_function = euclidean_distance,
                             engine = dense_engine()) {
  kernel <- structure(
    list(
      length_scale = wrapped_param(length_scale, length_scale_parametrisation),
      distance_function = resolve_distance_function(distance_function),
      engine = engine
    ),
    class = c("matern12_kernel", "stationary_kernel", "kernel")
  )
  validate_kernel(kernel)
  kernel
}

#' @export
pairwise.matern12_kernel <- function(kernel, x1, x2) {
  shape(kernel, kernel$distance_function(x1, x2))
}

#' @export
shape.matern12_kernel <- function(kernel, d) {
  length_scale <- unwrap_param(kernel$length_scale)
  exp(-d / length_scale)
}

#' @export
shape_grad.matern12_kernel <- function(kernel, d) {
  length_scale <- unwrap_param(kernel$length_scale)
  s <- shape(kernel, d)
  list(length_scale = s * d / length_scale^2)
}

#' @export
shape_grad_d.matern12_kernel <- function(kernel, d) {
  length_scale <- unwrap_param(kernel$length_scale)
  -shape(kernel, d) / length_scale
}

#' @rdname matern_kernels
#' @export
matern32_kernel <- function(length_scale,
                             length_scale_parametrisation = log_exp_parametrisation(),
                             distance_function = euclidean_distance,
                             engine = dense_engine()) {
  kernel <- structure(
    list(
      length_scale = wrapped_param(length_scale, length_scale_parametrisation),
      distance_function = resolve_distance_function(distance_function),
      engine = engine
    ),
    class = c("matern32_kernel", "stationary_kernel", "kernel")
  )
  validate_kernel(kernel)
  kernel
}

#' @export
pairwise.matern32_kernel <- function(kernel, x1, x2) {
  shape(kernel, kernel$distance_function(x1, x2))
}

#' @export
shape.matern32_kernel <- function(kernel, d) {
  length_scale <- unwrap_param(kernel$length_scale)
  scaled <- sqrt(3) * d / length_scale
  (1 + scaled) * exp(-scaled)
}

#' @export
shape_grad.matern32_kernel <- function(kernel, d) {
  length_scale <- unwrap_param(kernel$length_scale)
  scaled <- sqrt(3) * d / length_scale
  list(length_scale = scaled^2 * exp(-scaled) / length_scale)
}

#' @export
shape_grad_d.matern32_kernel <- function(kernel, d) {
  length_scale <- unwrap_param(kernel$length_scale)
  scaled <- sqrt(3) * d / length_scale
  -scaled * sqrt(3) * exp(-scaled) / length_scale
}

#' @rdname matern_kernels
#' @export
matern52_kernel <- function(length_scale,
                             length_scale_parametrisation = log_exp_parametrisation(),
                             distance_function = euclidean_distance,
                             engine = dense_engine()) {
  kernel <- structure(
    list(
      length_scale = wrapped_param(length_scale, length_scale_parametrisation),
      distance_function = resolve_distance_function(distance_function),
      engine = engine
    ),
    class = c("matern52_kernel", "stationary_kernel", "kernel")
  )
  validate_kernel(kernel)
  kernel
}

#' @export
pairwise.matern52_kernel <- function(kernel, x1, x2) {
  shape(kernel, kernel$distance_function(x1, x2))
}

#' @export
shape.matern52_kernel <- function(kernel, d) {
  length_scale <- unwrap_param(kernel$length_scale)
  scaled <- sqrt(5) * d / length_scale
  (1 + scaled + (5 / 3) * (d / length_scale)^2) * exp(-scaled)
}

#' @export
shape_grad.matern52_kernel <- function(kernel, d) {
  length_scale <- unwrap_param(kernel$length_scale)
  u <- d / length_scale
  scaled <- sqrt(5) * u
  list(length_scale = (5 / 3) * u^2 * (1 + scaled) * exp(-scaled) / length_scale)
}

#' @export
shape_grad_d.matern52_kernel <- function(kernel, d) {
  length_scale <- unwrap_param(kernel$length_scale)
  u <- d / length_scale
  scaled <- sqrt(5) * u
  -(5 / 3) * u * (1 + scaled) * exp(-scaled) / length_scale
}

#' Periodic kernel
#'
#' @description
#' \deqn{k(x_1, x_2) = \exp\left(-\frac{2 \sin^2(\pi d / \mathrm{period})}{\mathrm{length\_scale}^2}\right)}
#'
#' where \eqn{d} is `distance_function` (Euclidean by default). Both
#' `length_scale` and `period` are constrained to stay strictly positive.
#'
#' @section A very small `period`:
#' An extremely small `period` (e.g. `1e-300`) does not error -- `sin()`
#' still returns a finite, well-defined value -- but the result is numerical
#' aliasing: `d / period` is so large that floating-point rounding makes
#' `sin(pi * d / period)` behave essentially like noise as `d` varies, rather
#' than encoding any meaningful periodicity. The output stays finite and
#' symmetric, it is just not meaningful. No threshold is enforced in
#' `evaluate()` itself, since legitimate high-frequency use (a genuinely
#' tiny but intentional `period`) must not be blocked; see
#' `dev/ameliorations/README.md` (probe A10) for the measured behaviour.
#'
#' @param length_scale A strictly positive numeric value.
#' @param period A strictly positive numeric value.
#' @param length_scale_parametrisation A [parametrisation] for
#'   `length_scale`. Defaults to [log_exp_parametrisation()].
#' @param period_parametrisation A [parametrisation] for `period`. Defaults
#'   to [log_exp_parametrisation()].
#' @param distance_function A distance function, or a shortcut name (see
#'   [resolve_distance_function()]). Defaults to [euclidean_distance()].
#' @param engine A computation engine. Defaults to [dense_engine()].
#'
#' @return An object of class `periodic_kernel`/`stationary_kernel`/`kernel`.
#' @export
#' @examples
#' k <- periodic_kernel(length_scale = 1, period = 2)
#' evaluate(k, c(0, 1, 2))
periodic_kernel <- function(length_scale,
                             period,
                             length_scale_parametrisation = log_exp_parametrisation(),
                             period_parametrisation = log_exp_parametrisation(),
                             distance_function = euclidean_distance,
                             engine = dense_engine()) {
  kernel <- structure(
    list(
      length_scale = wrapped_param(length_scale, length_scale_parametrisation),
      period = wrapped_param(period, period_parametrisation),
      distance_function = resolve_distance_function(distance_function),
      engine = engine
    ),
    class = c("periodic_kernel", "stationary_kernel", "kernel")
  )
  validate_kernel(kernel)
  kernel
}

#' @export
pairwise.periodic_kernel <- function(kernel, x1, x2) {
  shape(kernel, kernel$distance_function(x1, x2))
}

#' @export
shape.periodic_kernel <- function(kernel, d) {
  length_scale <- unwrap_param(kernel$length_scale)
  period <- unwrap_param(kernel$period)
  exp(-2 * sin(pi * d / period)^2 / length_scale^2)
}

#' @export
shape_grad.periodic_kernel <- function(kernel, d) {
  length_scale <- unwrap_param(kernel$length_scale)
  period <- unwrap_param(kernel$period)
  arg <- pi * d / period
  s <- shape(kernel, d)
  list(
    length_scale = s * 4 * sin(arg)^2 / length_scale^3,
    period = s * 2 * pi * d * sin(2 * arg) / (length_scale^2 * period^2)
  )
}

#' @export
shape_grad_d.periodic_kernel <- function(kernel, d) {
  length_scale <- unwrap_param(kernel$length_scale)
  period <- unwrap_param(kernel$period)
  arg <- pi * d / period
  -2 * shape(kernel, d) * pi * sin(2 * arg) / (length_scale^2 * period)
}

#' Rational quadratic kernel
#'
#' @description
#' \deqn{k(x_1, x_2) = \left(1 + \frac{d}{2 \alpha \cdot \mathrm{length\_scale}^2}\right)^{-\alpha}}
#'
#' where \eqn{d} is `distance_function` (squared Euclidean by default).
#' `length_scale` and `alpha` are constrained to stay strictly positive.
#'
#' @param length_scale A strictly positive numeric value.
#' @param alpha A strictly positive numeric value.
#' @param length_scale_parametrisation A [parametrisation] for
#'   `length_scale`. Defaults to [log_exp_parametrisation()].
#' @param alpha_parametrisation A [parametrisation] for `alpha`. Defaults to
#'   [log_exp_parametrisation()].
#' @param distance_function A distance function, or a shortcut name (see
#'   [resolve_distance_function()]). Defaults to
#'   [squared_euclidean_distance()].
#' @param engine A computation engine. Defaults to [dense_engine()].
#'
#' @return An object of class `rational_quadratic_kernel`/
#'   `stationary_kernel`/`kernel`.
#' @export
#' @examples
#' k <- rational_quadratic_kernel(length_scale = 1, alpha = 1)
#' evaluate(k, c(0, 1, 2))
rational_quadratic_kernel <- function(length_scale,
                                       alpha,
                                       length_scale_parametrisation = log_exp_parametrisation(),
                                       alpha_parametrisation = log_exp_parametrisation(),
                                       distance_function = squared_euclidean_distance,
                                       engine = dense_engine()) {
  kernel <- structure(
    list(
      length_scale = wrapped_param(length_scale, length_scale_parametrisation),
      alpha = wrapped_param(alpha, alpha_parametrisation),
      distance_function = resolve_distance_function(distance_function),
      engine = engine
    ),
    class = c("rational_quadratic_kernel", "stationary_kernel", "kernel")
  )
  validate_kernel(kernel)
  kernel
}

#' @export
pairwise.rational_quadratic_kernel <- function(kernel, x1, x2) {
  shape(kernel, kernel$distance_function(x1, x2))
}

#' @export
shape.rational_quadratic_kernel <- function(kernel, d) {
  length_scale <- unwrap_param(kernel$length_scale)
  alpha <- unwrap_param(kernel$alpha)
  (1 + d / (2 * alpha * length_scale^2))^(-alpha)
}

#' @export
shape_grad.rational_quadratic_kernel <- function(kernel, d) {
  length_scale <- unwrap_param(kernel$length_scale)
  alpha <- unwrap_param(kernel$alpha)
  base <- 1 + d / (2 * alpha * length_scale^2)
  s <- base^(-alpha)
  list(
    length_scale = s * d / (length_scale^3 * base),
    alpha = s * ((base - 1) / base - log(base))
  )
}

#' @export
shape_grad_d.rational_quadratic_kernel <- function(kernel, d) {
  length_scale <- unwrap_param(kernel$length_scale)
  alpha <- unwrap_param(kernel$alpha)
  base <- 1 + d / (2 * alpha * length_scale^2)
  -shape(kernel, d) / (2 * length_scale^2 * base)
}

#' White noise kernel
#'
#' @description
#' Returns `noise` on the diagonal (where `x1 == x2`) and `0` elsewhere:
#'
#' \deqn{k(x_1, x_2) = \mathrm{noise} \cdot \mathbb{1}[x_1 = x_2]}
#'
#' `noise` is constrained to be *strictly* positive by the default
#' [log_exp_parametrisation()], which raises an error on a zero or negative
#' value rather than silently producing a degenerate kernel. Simply omit
#' `white_noise_kernel()` from a sum of kernels if no noise is desired.
#'
#' @param noise A strictly positive numeric value.
#' @param noise_parametrisation A [parametrisation] for `noise`. Defaults to
#'   [log_exp_parametrisation()].
#' @param engine A computation engine. Defaults to [dense_engine()].
#'
#' @return An object of class `white_noise_kernel`/`stationary_kernel`/
#'   `kernel`.
#' @export
#' @examples
#' k <- white_noise_kernel(noise = 0.1)
#' evaluate(k, c(0, 1, 2))
white_noise_kernel <- function(noise,
                                noise_parametrisation = log_exp_parametrisation(),
                                engine = dense_engine()) {
  kernel <- structure(
    list(
      noise = wrapped_param(noise, noise_parametrisation),
      distance_function = equality,
      engine = engine
    ),
    class = c("white_noise_kernel", "stationary_kernel", "kernel")
  )
  validate_kernel(kernel)
  kernel
}

#' @export
pairwise.white_noise_kernel <- function(kernel, x1, x2) {
  shape(kernel, kernel$distance_function(x1, x2))
}

#' @export
shape.white_noise_kernel <- function(kernel, d) {
  d * unwrap_param(kernel$noise)
}

#' @export
shape_grad.white_noise_kernel <- function(kernel, d) {
  list(noise = d)
}

#' @export
shape_grad_d.white_noise_kernel <- function(kernel, d) {
  0 * d + unwrap_param(kernel$noise)
}

#' Feature kernel
#'
#' @description
#' A kernel modelling the covariance between two latent "feature" dimensions
#' observed with their own length scales, plus a shared length scale `u`:
#'
#' \deqn{k(x_1, x_2) = \frac{\mathrm{variance}_1 \cdot \mathrm{variance}_2}
#' {(2\pi)^{D/2} \sqrt{\det \Sigma}}
#' \exp\left(-\frac{1}{2} \frac{d}{\Sigma}\right)}
#'
#' where \eqn{\Sigma = \mathrm{length\_scale}_1 + \mathrm{length\_scale}_2 +
#' \mathrm{length\_scale\_u}}, \eqn{D} is the input dimensionality, and
#' \eqn{d} is `distance_function` (squared Euclidean by default).
#'
#' `length_scale` and `variance` may each be given as a single value
#' (broadcast to both features) or a length-2 vector.
#'
#' @param length_scale A strictly positive numeric value or length-2 vector.
#' @param length_scale_u A strictly positive numeric value.
#' @param variance A strictly positive numeric value or length-2 vector.
#' @param length_scale_parametrisation,length_scale_u_parametrisation,variance_parametrisation
#'   [parametrisation] objects. Default to [log_exp_parametrisation()].
#' @param distance_function A distance function, or a shortcut name (see
#'   [resolve_distance_function()]). Defaults to
#'   [squared_euclidean_distance()].
#' @param engine A computation engine. Defaults to [dense_engine()].
#'
#' @return An object of class `feature_kernel`/`stationary_kernel`/`kernel`.
#' @export
#' @examples
#' k <- feature_kernel(length_scale = 1, length_scale_u = 1, variance = 1)
#' evaluate(k, c(0, 1, 2))
feature_kernel <- function(length_scale,
                            length_scale_u,
                            variance,
                            length_scale_parametrisation = log_exp_parametrisation(),
                            length_scale_u_parametrisation = log_exp_parametrisation(),
                            variance_parametrisation = log_exp_parametrisation(),
                            distance_function = squared_euclidean_distance,
                            engine = dense_engine()) {
  kernel <- structure(
    list(
      length_scale = wrapped_param(length_scale, length_scale_parametrisation),
      length_scale_u = wrapped_param(length_scale_u, length_scale_u_parametrisation),
      variance = wrapped_param(variance, variance_parametrisation),
      distance_function = resolve_distance_function(distance_function),
      engine = engine
    ),
    class = c("feature_kernel", "stationary_kernel", "kernel")
  )
  # Marked block_compatible now, ahead of block_kernel()'s own implementation
  # (build plan step 9), so the marker is not forgotten later -- see
  # dev/ameliorations/design/03-analyse-risques.Rmd, "Point de vigilance" on step 4.
  # feature_kernel() is explicitly designed to read a *pair* of values
  # (length_scale[1]/length_scale[2], variance[1]/variance[2]) rather than a
  # single scalar per hyperparameter, which is what block_kernel() requires.
  attr(kernel, "block_compatible") <- TRUE
  validate_kernel(kernel)
  kernel
}

#' @keywords internal
.feature_kernel_sigma <- function(kernel) {
  length_scale <- unwrap_param(kernel$length_scale)
  if (length(length_scale) == 1) length_scale <- rep(length_scale, 2)
  variance <- unwrap_param(kernel$variance)
  if (length(variance) == 1) variance <- rep(variance, 2)
  length_scale_u <- unwrap_param(kernel$length_scale_u)

  list(
    sigma_diag = length_scale[1] + length_scale[2] + length_scale_u,
    variance = variance
  )
}

#' @export
pairwise.feature_kernel <- function(kernel, x1, x2) {
  s <- .feature_kernel_sigma(kernel)
  sigma_det <- prod(s$sigma_diag)
  quadratic_form <- kernel$distance_function(x1, x2) / s$sigma_diag

  (s$variance[1] * s$variance[2] /
    ((2 * pi)^(length(x1) / 2) * sqrt(sigma_det))) *
    exp(-0.5 * quadratic_form)
}

#' @keywords internal
#' @exportS3Method
pairwise_matrix.feature_kernel <- function(kernel, X1, X2) {
  s <- .feature_kernel_sigma(kernel)
  sigma_det <- prod(s$sigma_diag)
  d <- distance_matrix(kernel$distance_function, X1, X2)

  (s$variance[1] * s$variance[2] /
    ((2 * pi)^(ncol(X1) / 2) * sqrt(sigma_det))) *
    exp(-0.5 * d / s$sigma_diag)
}

#' @keywords internal
#' @exportS3Method
pairwise_diag.feature_kernel <- function(kernel, X1, X2) {
  s <- .feature_kernel_sigma(kernel)
  sigma_det <- prod(s$sigma_diag)
  d <- diag_distance(kernel$distance_function, X1, X2)

  (s$variance[1] * s$variance[2] /
    ((2 * pi)^(ncol(X1) / 2) * sqrt(sigma_det))) *
    exp(-0.5 * d / s$sigma_diag)
}

# feature_kernel() has no shape() method (its formula isn't a function of
# the distance alone -- length_scale/length_scale_u also rescale the
# prefactor, not just the exponent), so its gradient is written directly
# against pairwise_matrix()/pairwise_diag() rather than via shape_grad().
# Writing k = A * Sigma^(-1/2) * exp(-0.5 d / Sigma) (A = variance[1] *
# variance[2] / (2 pi)^(D/2), constant w.r.t. Sigma) gives
# d(ln k)/d(Sigma) = -1/(2 Sigma) + 0.5 d / Sigma^2 = (d - Sigma) / (2
# Sigma^2); length_scale[1]/length_scale[2]/length_scale_u each contribute
# to Sigma with coefficient 1, so they share this same partial derivative.
# variance[1]/variance[2] enter k linearly, so d(k)/d(variance[i]) = k /
# variance[i].

#' @keywords internal
#' @exportS3Method
pairwise_matrix_grad.feature_kernel <- function(kernel, X1, X2) {
  s <- .feature_kernel_sigma(kernel)
  sigma <- s$sigma_diag
  k <- pairwise_matrix.feature_kernel(kernel, X1, X2)
  d <- distance_matrix(kernel$distance_function, X1, X2)
  grad_sigma <- k * (d - sigma) / (2 * sigma^2)
  list(
    length_scale = list(grad_sigma, grad_sigma),
    length_scale_u = grad_sigma,
    variance = list(k / s$variance[1], k / s$variance[2])
  )
}

#' @keywords internal
#' @exportS3Method
pairwise_diag_grad.feature_kernel <- function(kernel, X1, X2) {
  s <- .feature_kernel_sigma(kernel)
  sigma <- s$sigma_diag
  k <- pairwise_diag.feature_kernel(kernel, X1, X2)
  d <- diag_distance(kernel$distance_function, X1, X2)
  grad_sigma <- k * (d - sigma) / (2 * sigma^2)
  list(
    length_scale = list(grad_sigma, grad_sigma),
    length_scale_u = grad_sigma,
    variance = list(k / s$variance[1], k / s$variance[2])
  )
}
