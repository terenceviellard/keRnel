# kernel_grad() is checked two ways:
# - a couple of hand-derived closed-form values (independently verifiable,
#   same spirit as test-stationary-kernels.R's reference values)
# - broad coverage via central finite differences, the natural oracle for a
#   *derivative*: perturb_at() below bumps the idx-th trainable
#   hyperparameter (in kernel_grad_raw()'s own-fields-then-subkernels
#   order) by a small constrained-space delta, so it stays correct even
#   when two hyperparameters share a bare name (kupdate() would update
#   both at once, which a per-entry check must avoid).

.perturb_recurse <- function(kernel, idx, delta, component) {
  for (name in names(kernel)) {
    field <- kernel[[name]]
    if (inherits(field, "wrapped_param") && !isTRUE(field$frozen)) {
      idx <- idx - 1
      if (idx == 0) {
        val <- unwrap_param(field)
        if (!is.null(component)) val[component] <- val[component] + delta else val <- val + delta
        kernel[[name]]$wrapped <- wrap(field$parametrisation, val)
        return(list(kernel = kernel, idx = idx, found = TRUE))
      }
    }
  }
  if (!is.null(kernel$subkernels)) {
    for (i in seq_along(kernel$subkernels)) {
      sub_result <- .perturb_recurse(kernel$subkernels[[i]], idx, delta, component)
      idx <- sub_result$idx
      if (sub_result$found) {
        kernel$subkernels[[i]] <- sub_result$kernel
        return(list(kernel = kernel, idx = idx, found = TRUE))
      }
    }
  }
  list(kernel = kernel, idx = idx, found = FALSE)
}

perturb_at <- function(kernel, idx, delta, component = NULL) {
  .perturb_recurse(kernel, idx, delta, component)$kernel
}

expect_grad_matches_finite_diff <- function(kernel, x1, x2 = NULL, h = 1e-6, tolerance = 1e-4) {
  g <- kernel_grad(kernel, x1, x2)
  for (i in seq_along(g)) {
    entry <- g[[i]]
    comps <- if (is.list(entry)) seq_along(entry) else list(NULL)
    for (comp in comps) {
      analytic <- if (is.list(entry)) entry[[comp]] else entry
      k_plus <- perturb_at(kernel, i, h, component = comp)
      k_minus <- perturb_at(kernel, i, -h, component = comp)
      numeric_grad <- (evaluate(k_plus, x1, x2) - evaluate(k_minus, x1, x2)) / (2 * h)
      rel_err <- max(abs(numeric_grad - analytic)) / max(1, max(abs(analytic)))
      expect_lt(rel_err, tolerance)
    }
  }
}

x <- c(-1.3, 0, 0.7, 2.1)
x2 <- c(0.2, -0.5, 1.9, 3.0)

test_that("shape_grad.se_kernel matches the closed-form derivative at a known point", {
  # k(d) = exp(-0.5 d / l^2); dk/dl = k * d / l^3
  # squared_euclidean_distance(0, 1) == 1, l = 2
  k <- se_kernel(length_scale = 2)
  d <- 1
  expected <- exp(-0.5 * d / 4) * d / 8
  expect_equal(shape_grad(k, d)$length_scale, expected, tolerance = 1e-10)
})

test_that("kernel_grad() matches finite differences for every base kernel", {
  expect_grad_matches_finite_diff(se_kernel(length_scale = 1.7), x, x2)
  expect_grad_matches_finite_diff(matern12_kernel(length_scale = 0.9), x, x2)
  expect_grad_matches_finite_diff(matern32_kernel(length_scale = 1.3), x, x2)
  expect_grad_matches_finite_diff(matern52_kernel(length_scale = 2.2), x, x2)
  expect_grad_matches_finite_diff(periodic_kernel(length_scale = 1.1, period = 2.5), x, x2)
  expect_grad_matches_finite_diff(rational_quadratic_kernel(length_scale = 1.4, alpha = 0.8), x, x2)
  expect_grad_matches_finite_diff(white_noise_kernel(noise = 0.6), x, x)
  expect_grad_matches_finite_diff(linear_kernel(slope_var = 1.6), x, x2)
  expect_grad_matches_finite_diff(affine_kernel(slope_var = 1.6, offset = 0.5), x, x2)
  expect_grad_matches_finite_diff(polynomial_kernel(degree = 2, gamma = 0.7, constant = 1.2), x, x2)
  expect_grad_matches_finite_diff(polynomial_kernel(degree = 3, gamma = 0.7, constant = 1.2), x, x2)
  expect_grad_matches_finite_diff(sigmoid_kernel(alpha = 0.5, constant = 0.1), x, x2)
  expect_grad_matches_finite_diff(constant_kernel(value = 2.3), x, x2)
  expect_grad_matches_finite_diff(variance_kernel(variance = 1.9), x, x2)
})

test_that("kernel_grad() matches finite differences through sum_kernel()/product_kernel()", {
  expect_grad_matches_finite_diff(se_kernel(length_scale = 1.2) + white_noise_kernel(noise = 0.3), x, x2)
  expect_grad_matches_finite_diff(variance_kernel(variance = 2) * se_kernel(length_scale = 1.5), x, x2)
  expect_grad_matches_finite_diff(
    se_kernel(length_scale = 1.1) + variance_kernel(2) * matern32_kernel(length_scale = 0.8),
    x, x2
  )
})

test_that("kernel_grad() matches finite differences when a bare name collides across sub-kernels", {
  k <- se_kernel(length_scale = 1) + periodic_kernel(length_scale = 2, period = 3)
  expect_grad_matches_finite_diff(k, x, x2)
  expect_identical(sort(names(kernel_grad(k, x, x2))), sort(names(get_trainable_params(k))))
})

test_that("kernel_grad() matches finite differences under diagonal_engine()", {
  k <- se_kernel(length_scale = 1.4, engine = diagonal_engine())
  expect_grad_matches_finite_diff(k, x, x2)
  expect_true(is.numeric(kernel_grad(k, x, x2)$length_scale) && !is.matrix(kernel_grad(k, x, x2)$length_scale))
})

test_that("kernel_grad() matches finite differences for feature_kernel(), including vector hyperparameters", {
  X1 <- matrix(c(0, 1, 2, 0.5, -1, 1.5), ncol = 2)
  X2 <- matrix(c(0.3, 0.9, -0.4, 1.1, 2.0, -0.2), ncol = 2)
  k <- feature_kernel(length_scale = c(1.2, 0.9), length_scale_u = 0.5, variance = c(1.3, 0.8))
  expect_grad_matches_finite_diff(k, X1, X2)
  g <- kernel_grad(k, X1, X2)
  expect_length(g$length_scale, 2)
  expect_length(g$variance, 2)
})

test_that("kernel_grad() excludes frozen hyperparameters", {
  k <- freeze(se_kernel(length_scale = 1.5), "length_scale")
  expect_length(kernel_grad(k, x, x2), 0)
})

test_that("kernel_grad() names match get_trainable_params() names exactly", {
  k <- se_kernel(length_scale = 1) + periodic_kernel(length_scale = 2, period = 3) * white_noise_kernel(noise = 0.2)
  expect_identical(sort(names(kernel_grad(k, x, x2))), sort(names(get_trainable_params(k))))
})

test_that("kernel_grad() errors clearly on an unsupported (batching) wrapper kernel", {
  k <- batch_kernel(se_kernel(length_scale = c(1, 2)), batch_size = 2)
  expect_error(kernel_grad(k, x), "not implemented for `batch_kernel`")
})

test_that("kernel_grad() errors clearly when an unsupported kernel is composed with a supported one", {
  # Both sides produce a matching (n, n, 2) array, so evaluate() itself
  # succeeds -- the error must come from kernel_grad()'s own dispatch.
  k <- batch_kernel(se_kernel(length_scale = c(1, 2)), batch_size = 2) +
    batch_kernel(matern12_kernel(length_scale = c(1, 2)), batch_size = 2)
  expect_error(kernel_grad(k, x), "not implemented for `batch_kernel`")
})

test_that("kernel_grad() propagates evaluate()'s own validation (e.g. a shape mismatch)", {
  k <- se_kernel(length_scale = 1, engine = dense_engine()) +
    se_kernel(length_scale = 1, engine = diagonal_engine())
  expect_error(kernel_grad(k, x, x2), "dense_engine")
})

test_that("kernel_grad() matches finite differences through kernel_interactions()", {
  k1 <- se_kernel(length_scale = 1.3)
  k2 <- linear_kernel(slope_var = 0.9)
  k3 <- periodic_kernel(length_scale = 1.1, period = 2)
  expect_grad_matches_finite_diff(kernel_interactions(list(k1, k2, k3)), x)
})

test_that("kernel_grad() matches finite differences through exp_kernel()/log_kernel()", {
  expect_grad_matches_finite_diff(exp_kernel(se_kernel(length_scale = 1.2)), x, x2)
  expect_grad_matches_finite_diff(log_kernel(variance_kernel(variance = 2.5) * se_kernel(length_scale = 1.1)), x, x2)
})

test_that("kernel_grad() matches finite differences through active_dims_kernel()", {
  X1 <- matrix(c(0, 10, -1, 1, 11, 2, 2, 12, -3, -1.5, 9, 0.5), ncol = 3, byrow = TRUE)
  X2 <- matrix(c(0.3, 9.5, 1, -0.7, 10.2, -1, 1.8, 11.5, 2.5), ncol = 3, byrow = TRUE)
  k <- active_dims_kernel(se_kernel(length_scale = 1.4) + white_noise_kernel(noise = 0.2), active_dims = c(2, 3))
  expect_grad_matches_finite_diff(k, X1, X2)
})

test_that("kernel_grad() matches finite differences through ard_kernel(), across supported distance_functions", {
  X1 <- matrix(c(0, 1, 2, 3, 0.5, -1, 1.5, 2.5, -0.5, 0.3, 1.1, -0.8), ncol = 3, byrow = TRUE)
  X2 <- matrix(c(0.4, -0.2, 1.3, 2.1, 0.9, -1.4, 0.2, 1.6, -0.3, 1.2, 0.7, -0.5), ncol = 3, byrow = TRUE)

  # squared_euclidean_distance() (se_kernel's default)
  expect_grad_matches_finite_diff(ard_kernel(se_kernel, length_scales = c(0.8, 1.5, 2.2)), X1, X2)

  # euclidean_distance() (matern32_kernel's default), plus an extra free
  # inner hyperparameter alongside the frozen ARD one
  expect_grad_matches_finite_diff(
    ard_kernel(periodic_kernel, length_scales = c(1.2, 0.6), period = 3),
    X1[, 1:2], X2[, 1:2]
  )

  # manhattan_distance()
  expect_grad_matches_finite_diff(
    ard_kernel(matern12_kernel, length_scales = c(1.1, 0.8), distance_function = "manhattan"),
    X1[, 1:2], X2[, 1:2]
  )

  # dot_product() (linear_kernel's default)
  expect_grad_matches_finite_diff(
    ard_kernel(linear_kernel, length_scales = c(1.3, 0.9, 2.1), active = "slope_var"),
    X1, X2
  )
})

test_that("kernel_grad() matches finite differences with ard_kernel() nested inside other wrappers/operators", {
  X1 <- matrix(c(0, 1, 2, 3, 0.5, -1, 1.5, 2.5, -0.5, 0.3, 1.1, -0.8), ncol = 3, byrow = TRUE)
  X2 <- matrix(c(0.4, -0.2, 1.3, 2.1, 0.9, -1.4, 0.2, 1.6, -0.3, 1.2, 0.7, -0.5), ncol = 3, byrow = TRUE)

  expect_grad_matches_finite_diff(
    ard_kernel(se_kernel, length_scales = c(1, 1.5)) + white_noise_kernel(noise = 0.3),
    X1[, 1:2], X2[, 1:2]
  )
  expect_grad_matches_finite_diff(
    active_dims_kernel(ard_kernel(se_kernel, length_scales = c(1, 2)), active_dims = c(1, 3)),
    X1, X2
  )
})

test_that("ard_kernel() gradient errors clearly for an unsupported distance_function", {
  k <- ard_kernel(se_kernel, length_scales = c(1, 1), distance_function = "cosine")
  X <- matrix(c(0, 1, 1, 0, 1, 1, 2, 0), ncol = 2, byrow = TRUE)
  expect_error(kernel_grad(k, X), "no analytic gradient for this `distance_function`")
})

test_that("ard_kernel() gradient errors clearly when the wrapped kernel has no shape()", {
  k <- ard_kernel(feature_kernel, length_scales = c(1, 1),
                   active = "length_scale_u", length_scale = c(1, 1), variance = c(1, 1))
  X <- matrix(c(0, 1, 1, 0, 1, 1, 2, 0), ncol = 2, byrow = TRUE)
  expect_error(kernel_grad(k, X), "not implemented for ard_kernel\\(feature_kernel")
})

# kernel_grad_free() is checked positionally against a central finite
# difference taken directly in free (wrapped) space, via
# get_free_params()/set_free_params() -- the space it actually promises to
# match, independent of kernel_grad()/kupdate() and any name collision.
expect_free_grad_matches_finite_diff <- function(kernel, x1, x2 = NULL, h = 1e-6, tolerance = 1e-4) {
  g <- kernel_grad_free(kernel, x1, x2)
  free0 <- get_free_params(kernel)
  expect_length(g, length(free0))
  for (i in seq_along(g)) {
    f_plus <- free0
    f_plus[i] <- f_plus[i] + h
    f_minus <- free0
    f_minus[i] <- f_minus[i] - h
    numeric_grad <- (evaluate(set_free_params(kernel, f_plus), x1, x2) -
      evaluate(set_free_params(kernel, f_minus), x1, x2)) / (2 * h)
    rel_err <- max(abs(numeric_grad - g[[i]])) / max(1, max(abs(g[[i]])))
    expect_lt(rel_err, tolerance)
  }
}

test_that("kernel_grad_free() matches finite differences in free space, including a name collision", {
  expect_free_grad_matches_finite_diff(se_kernel(length_scale = 1.7), x, x2)
  expect_free_grad_matches_finite_diff(se_kernel(length_scale = 1) + se_kernel(length_scale = 3), x, x2)
  expect_free_grad_matches_finite_diff(variance_kernel(variance = 2) * se_kernel(length_scale = 1.5), x, x2)
  expect_free_grad_matches_finite_diff(constant_kernel(value = -1.5), x, x2)
})

test_that("kernel_grad_free() matches finite differences under a non-default parametrisation", {
  expect_free_grad_matches_finite_diff(
    se_kernel(length_scale = 3, length_scale_parametrisation = bounded_parametrisation(0, 10)), x, x2
  )
  expect_free_grad_matches_finite_diff(
    se_kernel(length_scale = 1.4, length_scale_parametrisation = softplus_parametrisation()), x, x2
  )
  expect_free_grad_matches_finite_diff(
    se_kernel(
      length_scale = 2,
      length_scale_parametrisation = parametrisation_chain(bounded_parametrisation(0, 10))
    ),
    x, x2
  )
})

test_that("kernel_grad_free() matches finite differences for vector-valued hyperparameters", {
  X1 <- matrix(c(0, 1, 2, 0.5, -1, 1.5), ncol = 2)
  X2 <- matrix(c(0.3, 0.9, -0.4, 1.1, 2.0, -0.2), ncol = 2)
  fk <- feature_kernel(length_scale = c(1.2, 0.9), length_scale_u = 0.5, variance = c(1.3, 0.8))
  expect_free_grad_matches_finite_diff(fk, X1, X2)

  Xa <- matrix(c(0, 1, 2, 3, 0.5, -1, 1.5, 2.5, -0.5, 0.3, 1.1, -0.8), ncol = 3, byrow = TRUE)
  Xb <- matrix(c(0.4, -0.2, 1.3, 2.1, 0.9, -1.4, 0.2, 1.6, -0.3, 1.2, 0.7, -0.5), ncol = 3, byrow = TRUE)
  expect_free_grad_matches_finite_diff(ard_kernel(se_kernel, length_scales = c(0.8, 1.5, 2.2)), Xa, Xb)
})

test_that("kernel_grad_free() excludes frozen hyperparameters, in step with get_free_params()", {
  k <- freeze(se_kernel(length_scale = 1.5) + white_noise_kernel(noise = 0.2), "length_scale")
  g <- kernel_grad_free(k, x, x2)
  expect_length(g, length(get_free_params(k)))
  expect_length(g, 1)
})
