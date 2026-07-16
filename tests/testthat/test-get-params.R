test_that("get_params() returns the kernel's hyperparameters by name", {
  k <- se_kernel(length_scale = 2.5)
  params <- get_params(k)
  expect_equal(params[["length_scale"]], 2.5, tolerance = 1e-8)
})

test_that("get_trainable_params() excludes frozen parameters", {
  k <- se_kernel(length_scale = 2.5)
  k_frozen <- freeze(k, "length_scale")

  expect_equal(get_params(k_frozen)[["length_scale"]], 2.5, tolerance = 1e-8)
  expect_length(get_trainable_params(k_frozen), 0)
})

test_that("flatten_params() returns a list of wrapped_param objects", {
  k <- se_kernel(length_scale = 1)
  flat <- flatten_params(k)
  expect_named(flat, "length_scale")
  expect_true(inherits(flat$length_scale, "wrapped_param"))
})

# Name disambiguation (finding #6): get_params()/get_trainable_params() must
# never silently return a vector with duplicate names.

test_that("get_params() leaves names untouched when there is no collision", {
  k <- se_kernel(length_scale = 1) + linear_kernel(slope_var = 2)
  params <- get_params(k)
  expect_named(params, c("length_scale", "slope_var"), ignore.order = TRUE)
})

test_that("get_params() qualifies a shared name by owning class", {
  k <- se_kernel(length_scale = 1) + periodic_kernel(length_scale = 2, period = 3)
  params <- get_params(k)
  expect_setequal(
    names(params),
    c("length_scale (se_kernel)", "length_scale (periodic_kernel)", "period")
  )
  expect_equal(params[["length_scale (se_kernel)"]], 1, tolerance = 1e-8)
  expect_equal(params[["length_scale (periodic_kernel)"]], 2, tolerance = 1e-8)
})

test_that("get_params() numbers occurrences when the same class repeats", {
  k <- se_kernel(length_scale = 1) + se_kernel(length_scale = 2)
  params <- get_params(k)
  expect_setequal(names(params), c("length_scale (se_kernel #1)", "length_scale (se_kernel #2)"))
  expect_equal(unname(sort(params)), c(1, 2), tolerance = 1e-8)
})

test_that("get_trainable_params() disambiguates the same way as get_params()", {
  k <- se_kernel(length_scale = 1) + se_kernel(length_scale = 2)
  expect_setequal(
    names(get_trainable_params(k)),
    c("length_scale (se_kernel #1)", "length_scale (se_kernel #2)")
  )
})

test_that("get_params() disambiguation does not disturb kupdate()/freeze() (bare-name propagation)", {
  k <- se_kernel(length_scale = 1) + se_kernel(length_scale = 2)
  k2 <- kupdate(k, length_scale = 9)
  expect_equal(unname(get_params(k2)), c(9, 9), tolerance = 1e-8)
})

# get_free_params()/set_free_params() (item 1b): the unconstrained (wrapped)
# space, positional rather than name-keyed.

test_that("get_free_params() returns the wrapped (unconstrained) value, not the constrained one", {
  k <- se_kernel(length_scale = 2)
  expect_equal(get_free_params(k), log(2), tolerance = 1e-12) # log_exp_parametrisation
  expect_equal(get_trainable_params(k)[["length_scale"]], 2, tolerance = 1e-12)
})

test_that("get_free_params() excludes frozen parameters, like get_trainable_params()", {
  k <- freeze(se_kernel(length_scale = 2), "length_scale")
  expect_length(get_free_params(k), 0)
})

test_that("get_free_params() is unnamed and positional", {
  k <- se_kernel(length_scale = 1) + se_kernel(length_scale = 2)
  expect_null(names(get_free_params(k)))
  expect_length(get_free_params(k), 2)
})

test_that("set_free_params()/get_free_params() round-trip bit-exactly", {
  k <- se_kernel(length_scale = 1) + white_noise_kernel(noise = 0.1)
  free <- get_free_params(k)
  k2 <- set_free_params(k, free)
  expect_identical(get_free_params(k2), free)
  expect_identical(k2, k)
})

test_that("set_free_params() can set two same-named hyperparameters independently (unlike kupdate())", {
  k <- se_kernel(length_scale = 1) + se_kernel(length_scale = 2)
  k2 <- set_free_params(k, c(log(9), log(3)))
  expect_equal(unname(get_params(k2)), c(9, 3), tolerance = 1e-8)
})

test_that("set_free_params() leaves frozen parameters untouched", {
  k <- freeze(se_kernel(length_scale = 1) + white_noise_kernel(noise = 0.1), "length_scale")
  k2 <- set_free_params(k, get_free_params(k) + 1)
  expect_equal(get_params(k2)[["length_scale"]], 1, tolerance = 1e-8) # frozen, unchanged
  expect_false(isTRUE(all.equal(get_params(k2)[["noise"]], 0.1, tolerance = 1e-8)))
})

test_that("set_free_params() handles a multi-valued hyperparameter as one contiguous block", {
  k <- ard_kernel(se_kernel, length_scales = c(1, 2, 3))
  free <- get_free_params(k)
  k2 <- set_free_params(k, free + 1)
  expect_equal(get_params(k2)[["length_scales1"]], exp(log(1) + 1), tolerance = 1e-8)
  expect_equal(get_params(k2)[["length_scales3"]], exp(log(3) + 1), tolerance = 1e-8)
})

test_that("set_free_params() errors on the wrong length", {
  k <- se_kernel(length_scale = 1)
  expect_error(set_free_params(k, c(1, 2)), "length 1")
})

test_that("set_free_params() errors on non-finite values", {
  k <- se_kernel(length_scale = 1)
  expect_error(set_free_params(k, Inf), "finite")
  expect_error(set_free_params(k, NA_real_), "finite")
})
