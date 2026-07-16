# Reference values below are derived by hand from the SE formula
# k(x1, x2) = exp(-0.5 * squared_euclidean_distance(x1, x2) / length_scale^2)
# (see dev/ameliorations/design/03-analyse-risques.Rmd, risk #5: tests must be independently
# verifiable, since we deliberately do not rely on externally generated
# fixtures).

test_that("se_kernel equals 1 at zero distance, for any length_scale", {
  k1 <- se_kernel(length_scale = 1)
  k2 <- se_kernel(length_scale = 7.3)
  expect_equal(pairwise(k1, c(0, 0), c(0, 0)), 1)
  expect_equal(pairwise(k2, c(3, -1), c(3, -1)), 1)
})

test_that("se_kernel matches the closed-form value for a known distance", {
  # squared_euclidean_distance(c(0,0), c(1,0)) == 1, length_scale = 1
  # k = exp(-0.5 * 1 / 1^2) = exp(-0.5)
  k <- se_kernel(length_scale = 1)
  expect_equal(pairwise(k, c(0, 0), c(1, 0)), exp(-0.5), tolerance = 1e-8)
})

test_that("se_kernel is symmetric", {
  k <- se_kernel(length_scale = 1.7)
  expect_equal(pairwise(k, c(1, 2), c(3, -1)), pairwise(k, c(3, -1), c(1, 2)))
})

test_that("se_kernel decreases monotonically with squared distance", {
  k <- se_kernel(length_scale = 1)
  near <- pairwise(k, c(0, 0), c(0.5, 0))
  far <- pairwise(k, c(0, 0), c(2, 0))
  expect_gt(near, far)
})

test_that("se_kernel's full evaluate() produces a symmetric positive semi-definite Gram matrix", {
  k <- se_kernel(length_scale = 1)
  x <- c(-1, 0, 0.5, 2)
  K <- evaluate(k, x)

  expect_equal(K, t(K), tolerance = 1e-10)
  eigenvalues <- eigen(K, symmetric = TRUE, only.values = TRUE)$values
  expect_true(all(eigenvalues > -1e-8))
})

test_that("se_kernel handles multi-dimensional inputs", {
  k <- se_kernel(length_scale = 1)
  X <- matrix(c(0, 0, 1, 0, 0, 1), nrow = 3, byrow = TRUE)
  K <- evaluate(k, X)
  expect_equal(dim(K), c(3, 3))
  expect_equal(diag(K), c(1, 1, 1))
})

test_that("se_kernel rejects non-positive length_scale at construction", {
  expect_error(se_kernel(length_scale = 0), "strictly positive")
  expect_error(se_kernel(length_scale = -1), "strictly positive")
})

test_that("se_kernel accepts a distance_function shortcut name", {
  k <- se_kernel(length_scale = 1, distance_function = "euclidean")
  expect_identical(k$distance_function, euclidean_distance)
})

test_that("rbf_kernel is an alias for se_kernel", {
  expect_identical(rbf_kernel, se_kernel)
})

test_that("se_kernel's class vector reflects its category hierarchy", {
  k <- se_kernel(length_scale = 1)
  expect_equal(class(k), c("se_kernel", "stationary_kernel", "kernel"))
})

# matern12_kernel: k(x1, x2) = exp(-euclidean_distance(x1, x2) / length_scale)

test_that("matern12_kernel equals 1 at zero distance", {
  k <- matern12_kernel(length_scale = 2)
  expect_equal(pairwise(k, c(1, 1), c(1, 1)), 1)
})

test_that("matern12_kernel matches the closed-form value for a known distance", {
  # euclidean_distance(c(0,0), c(3,4)) == 5, length_scale = 2
  # k = exp(-5 / 2)
  k <- matern12_kernel(length_scale = 2)
  expect_equal(pairwise(k, c(0, 0), c(3, 4)), exp(-2.5), tolerance = 1e-8)
})

test_that("matern12_kernel's evaluate() agrees with the vectorised fast path", {
  k <- matern12_kernel(length_scale = 1.5)
  x <- c(-1, 0, 0.5, 2)
  expect_equal(
    evaluate(k, x),
    map_pairs(function(a, b) pairwise(k, a, b), as_matrix_input(x), as_matrix_input(x)),
    tolerance = 1e-10
  )
})

# matern32_kernel: k = (1 + sqrt(3)*r/l) * exp(-sqrt(3)*r/l), r = euclidean_distance

test_that("matern32_kernel equals 1 at zero distance", {
  k <- matern32_kernel(length_scale = 2)
  expect_equal(pairwise(k, c(1, 1), c(1, 1)), 1)
})

test_that("matern32_kernel matches the closed-form value for a known distance", {
  # euclidean_distance(c(0,0), c(3,4)) == 5, length_scale = 2
  # scaled = sqrt(3) * 5 / 2; k = (1 + scaled) * exp(-scaled)
  k <- matern32_kernel(length_scale = 2)
  scaled <- sqrt(3) * 5 / 2
  expect_equal(pairwise(k, c(0, 0), c(3, 4)), (1 + scaled) * exp(-scaled), tolerance = 1e-8)
})

test_that("matern32_kernel's evaluate() agrees with the vectorised fast path", {
  k <- matern32_kernel(length_scale = 0.8)
  x <- c(-1, 0, 0.5, 2)
  expect_equal(
    evaluate(k, x),
    map_pairs(function(a, b) pairwise(k, a, b), as_matrix_input(x), as_matrix_input(x)),
    tolerance = 1e-10
  )
})

# matern52_kernel: k = (1 + sqrt(5)*r/l + (5/3)*(r/l)^2) * exp(-sqrt(5)*r/l)

test_that("matern52_kernel equals 1 at zero distance", {
  k <- matern52_kernel(length_scale = 2)
  expect_equal(pairwise(k, c(1, 1), c(1, 1)), 1)
})

test_that("matern52_kernel matches the closed-form value for a known distance", {
  # euclidean_distance(c(0,0), c(3,4)) == 5, length_scale = 2
  k <- matern52_kernel(length_scale = 2)
  scaled <- sqrt(5) * 5 / 2
  expected <- (1 + scaled + (5 / 3) * (5 / 2)^2) * exp(-scaled)
  expect_equal(pairwise(k, c(0, 0), c(3, 4)), expected, tolerance = 1e-8)
})

test_that("matern52_kernel's evaluate() agrees with the vectorised fast path", {
  k <- matern52_kernel(length_scale = 1.1)
  x <- c(-1, 0, 0.5, 2)
  expect_equal(
    evaluate(k, x),
    map_pairs(function(a, b) pairwise(k, a, b), as_matrix_input(x), as_matrix_input(x)),
    tolerance = 1e-10
  )
})

test_that("Matern kernels rank smoothness-induced decay consistently with length_scale", {
  # All three Matern kernels must equal 1 at zero distance and decay strictly
  # as distance increases, for any smoothness.
  for (k in list(matern12_kernel(1), matern32_kernel(1), matern52_kernel(1))) {
    near <- pairwise(k, c(0, 0), c(0.5, 0))
    far <- pairwise(k, c(0, 0), c(2, 0))
    expect_gt(near, far)
  }
})

# periodic_kernel: k = exp(-2 * sin(pi * d / period)^2 / length_scale^2), d = euclidean_distance

test_that("periodic_kernel equals 1 at zero distance", {
  k <- periodic_kernel(length_scale = 1, period = 2)
  expect_equal(pairwise(k, c(0, 0), c(0, 0)), 1)
})

test_that("periodic_kernel equals 1 at integer multiples of the period", {
  k <- periodic_kernel(length_scale = 1, period = 2)
  expect_equal(pairwise(k, c(0), c(2)), 1, tolerance = 1e-10)
  expect_equal(pairwise(k, c(0), c(4)), 1, tolerance = 1e-10)
})

test_that("periodic_kernel matches the closed-form value for a known distance", {
  # d = 1, period = 2, length_scale = 1
  # k = exp(-2 * sin(pi/2)^2 / 1) = exp(-2)
  k <- periodic_kernel(length_scale = 1, period = 2)
  expect_equal(pairwise(k, c(0), c(1)), exp(-2), tolerance = 1e-8)
})

test_that("periodic_kernel's evaluate() agrees with the vectorised fast path", {
  k <- periodic_kernel(length_scale = 1.3, period = 3)
  x <- c(-1, 0, 0.5, 2)
  expect_equal(
    evaluate(k, x),
    map_pairs(function(a, b) pairwise(k, a, b), as_matrix_input(x), as_matrix_input(x)),
    tolerance = 1e-10
  )
})

test_that("periodic_kernel rejects non-positive period", {
  expect_error(periodic_kernel(length_scale = 1, period = 0), "strictly positive")
  expect_error(periodic_kernel(length_scale = 1, period = -1), "strictly positive")
})

# rational_quadratic_kernel: k = (1 + d / (2*alpha*l^2))^(-alpha), d = squared_euclidean_distance

test_that("rational_quadratic_kernel equals 1 at zero distance", {
  k <- rational_quadratic_kernel(length_scale = 1, alpha = 2)
  expect_equal(pairwise(k, c(0, 0), c(0, 0)), 1)
})

test_that("rational_quadratic_kernel matches the closed-form value for a known distance", {
  # squared_euclidean_distance(c(0,0), c(1,0)) == 1, length_scale = 1, alpha = 2
  # k = (1 + 1 / (2*2*1))^(-2) = (1.25)^(-2)
  k <- rational_quadratic_kernel(length_scale = 1, alpha = 2)
  expect_equal(pairwise(k, c(0, 0), c(1, 0)), 1.25^(-2), tolerance = 1e-8)
})

test_that("rational_quadratic_kernel's evaluate() agrees with the vectorised fast path", {
  k <- rational_quadratic_kernel(length_scale = 0.9, alpha = 1.5)
  x <- c(-1, 0, 0.5, 2)
  expect_equal(
    evaluate(k, x),
    map_pairs(function(a, b) pairwise(k, a, b), as_matrix_input(x), as_matrix_input(x)),
    tolerance = 1e-10
  )
})

test_that("rational_quadratic_kernel approaches se_kernel-like decay as alpha grows", {
  # As alpha -> Inf, RQ(l, alpha) -> SE(l). Sanity check at a large alpha.
  rq <- rational_quadratic_kernel(length_scale = 1, alpha = 1e6)
  se <- se_kernel(length_scale = 1)
  expect_equal(pairwise(rq, c(0, 0), c(0.3, 0.1)), pairwise(se, c(0, 0), c(0.3, 0.1)), tolerance = 1e-4)
})

# white_noise_kernel: k = noise * equality(x1, x2)

test_that("white_noise_kernel returns noise on the diagonal and 0 off it", {
  k <- white_noise_kernel(noise = 0.5)
  expect_equal(pairwise(k, c(1, 2), c(1, 2)), 0.5)
  expect_equal(pairwise(k, c(1, 2), c(1, 3)), 0)
})

test_that("white_noise_kernel's evaluate() produces a diagonal matrix", {
  k <- white_noise_kernel(noise = 2)
  x <- c(1, 2, 3)
  K <- evaluate(k, x)
  expect_equal(diag(K), c(2, 2, 2))
  expect_equal(K[upper.tri(K)], rep(0, 3))
})

test_that("white_noise_kernel rejects non-positive noise", {
  expect_error(white_noise_kernel(noise = 0), "strictly positive")
  expect_error(white_noise_kernel(noise = -1), "strictly positive")
})

# feature_kernel: a scaled SE-like kernel, see R/stationary-kernels.R

test_that("feature_kernel matches the closed-form value for a known distance", {
  # length_scale = c(1, 1), length_scale_u = 1, variance = c(2, 3)
  # sigma_diag = 1 + 1 + 1 = 3, sigma_det = 3
  # squared_euclidean_distance(c(0,0), c(1,0)) = 1
  # k = (2*3) / ((2*pi)^1 * sqrt(3)) * exp(-0.5 * 1/3)
  k <- feature_kernel(length_scale = 1, length_scale_u = 1, variance = c(2, 3))
  expected <- (2 * 3) / ((2 * pi)^1 * sqrt(3)) * exp(-0.5 * 1 / 3)
  expect_equal(pairwise(k, c(0, 0), c(1, 0)), expected, tolerance = 1e-8)
})

test_that("feature_kernel is marked block_compatible ahead of block_kernel()'s implementation", {
  # See dev/ameliorations/design/03-analyse-risques.Rmd, "Point de vigilance" on step 4: this
  # marker must be set now so it isn't forgotten by the time block_kernel()
  # is implemented (build plan step 9).
  k <- feature_kernel(length_scale = 1, length_scale_u = 1, variance = 1)
  expect_true(isTRUE(attr(k, "block_compatible")))
})

test_that("feature_kernel's evaluate() agrees with the vectorised fast path", {
  k <- feature_kernel(length_scale = c(0.7, 1.3), length_scale_u = 0.5, variance = c(1.5, 2))
  X <- matrix(c(0, 0, 1, 0, 0, 1, -1, 2), ncol = 2, byrow = TRUE)
  expect_equal(
    evaluate(k, X),
    map_pairs(function(a, b) pairwise(k, a, b), X, X),
    tolerance = 1e-10
  )
})
