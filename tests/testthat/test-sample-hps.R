test_that("sample_hps() resamples a hyperparameter within the bounds of a uniform prior", {
  set.seed(1)
  k <- se_kernel(length_scale = 1)
  k2 <- sample_hps(k, priors = list(length_scale = c(0.1, 5)))
  value <- unwrap_param(k2$length_scale)
  expect_true(value >= 0.1 && value <= 5)
})

test_that("sample_hps() changes the value (with overwhelming probability)", {
  set.seed(1)
  k <- se_kernel(length_scale = 1)
  k2 <- sample_hps(k, priors = list(length_scale = c(0.1, 5)))
  expect_false(isTRUE(all.equal(unwrap_param(k$length_scale), unwrap_param(k2$length_scale))))
})

test_that("sample_hps() supports normal priors via distribution = stats::rnorm", {
  set.seed(42)
  k <- linear_kernel(slope_var = 1)
  k2 <- sample_hps(k, priors = list(slope_var = c(2, 1e-6)), distribution = stats::rnorm)
  expect_equal(unwrap_param(k2$slope_var), 2, tolerance = 1e-4)
})

test_that("sample_hps() leaves unspecified hyperparameters untouched", {
  set.seed(1)
  k <- periodic_kernel(length_scale = 1, period = 7)
  k2 <- sample_hps(k, priors = list(length_scale = c(0.1, 5)))
  expect_equal(unwrap_param(k2$period), 7)
})

test_that("sample_hps() skips a frozen hyperparameter", {
  set.seed(1)
  k <- freeze(se_kernel(length_scale = 1), "length_scale")
  k2 <- sample_hps(k, priors = list(length_scale = c(0.1, 5)))
  expect_equal(unwrap_param(k2$length_scale), 1)
  expect_true(k2$length_scale$frozen)
})

test_that("sample_hps() recurses through $subkernels of a composite kernel", {
  set.seed(1)
  k <- se_kernel(length_scale = 1) + white_noise_kernel(noise = 1)
  k2 <- sample_hps(k, priors = list(length_scale = c(0.1, 5), noise = c(0.01, 0.5)))
  ls <- unwrap_param(k2$subkernels[[1]]$length_scale)
  noise <- unwrap_param(k2$subkernels[[2]]$noise)
  expect_true(ls >= 0.1 && ls <= 5)
  expect_true(noise >= 0.01 && noise <= 0.5)
})

test_that("sample_hps() resamples each element of a vector-valued hyperparameter independently", {
  set.seed(1)
  k <- ard_kernel(se_kernel, length_scales = c(1, 1, 1))
  k2 <- sample_hps(k, priors = list(length_scales = c(0.1, 10)))
  values <- unwrap_param(k2$length_scales)
  expect_length(values, 3)
  expect_true(all(values >= 0.1 & values <= 10))
})

test_that("sample_hps() errors on a prior name that matches no hyperparameter", {
  k <- se_kernel(length_scale = 1)
  expect_error(sample_hps(k, priors = list(not_a_param = c(0, 1))), "No such hyperparameter")
})

test_that("sample_hps() propagates a constraint violation from wrap()", {
  set.seed(7)
  k <- se_kernel(length_scale = 1)
  # A uniform prior straddling 0 will, eventually, draw a non-positive value,
  # which log_exp_parametrisation()'s wrap() must reject.
  expect_error(
    sample_hps(k, priors = list(length_scale = c(-1, -0.5))),
    "strictly positive"
  )
})

test_that("sample_hps() with an empty priors list returns the kernel unchanged", {
  k <- se_kernel(length_scale = 3)
  k2 <- sample_hps(k, priors = list())
  expect_identical(k, k2)
})
