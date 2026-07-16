test_that("kupdate() updates the named hyperparameter", {
  k <- se_kernel(length_scale = 1)
  k2 <- kupdate(k, length_scale = 5)
  expect_equal(get_params(k2)[["length_scale"]], 5, tolerance = 1e-8)
})

test_that("kupdate() does not mutate the original kernel (copy-on-modify)", {
  k <- se_kernel(length_scale = 1)
  kupdate(k, length_scale = 5)
  expect_equal(get_params(k)[["length_scale"]], 1, tolerance = 1e-8)
})

test_that("kupdate() re-validates the new value via the parametrisation's wrap()", {
  k <- se_kernel(length_scale = 1)
  expect_error(kupdate(k, length_scale = -1), "strictly positive")
})

test_that("kupdate() errors on an unknown hyperparameter name (typo protection)", {
  k <- se_kernel(length_scale = 1)
  expect_error(kupdate(k, lenght_scale = 5), "No such hyperparameter")
})

test_that("kupdate() with no arguments returns the kernel unchanged", {
  k <- se_kernel(length_scale = 1)
  expect_equal(kupdate(k), k)
})

test_that("kupdate() errors rather than silently overwriting a frozen hyperparameter", {
  k <- freeze(se_kernel(length_scale = 1), "length_scale")
  expect_error(kupdate(k, length_scale = 5), "frozen")
  expect_equal(get_params(k)[["length_scale"]], 1, tolerance = 1e-8)
})

test_that("kupdate() respects ard_kernel()'s frozen inner length_scale (finding #3)", {
  k <- ard_kernel(se_kernel, length_scales = c(1, 2, 3))
  expect_error(kupdate(k, length_scale = 99), "frozen")
  # The protected inner length_scale is untouched; length_scales is a
  # different name and remains updatable normally.
  k2 <- kupdate(k, length_scales = c(4, 5, 6))
  expect_equal(unname(get_params(k2)), c(4, 5, 6, 1), tolerance = 1e-10)
})

test_that("kupdate() updates only the non-frozen occurrences when a name is shared", {
  frozen_part <- freeze(se_kernel(length_scale = 1), "length_scale")
  trainable_part <- matern32_kernel(length_scale = 2)
  k <- sum_kernel(frozen_part, trainable_part)

  k2 <- kupdate(k, length_scale = 9)
  vals <- unname(get_params(k2))
  expect_equal(sort(vals), sort(c(1, 9)), tolerance = 1e-8)
})

test_that("kupdate() refuses to silently collapse a multi-valued (batch-axis) hyperparameter to a shared scalar", {
  bk <- batch_kernel(se_kernel(length_scale = c(1, 2, 3)), batch_size = 3)
  expect_error(kupdate(bk, length_scale = 99), "3 value")
  # the kernel is genuinely untouched -- no partial update happened
  expect_equal(unname(get_params(bk)[grepl("length_scale", names(get_params(bk)))]), c(1, 2, 3), tolerance = 1e-8)
})

test_that("kupdate() still allows replacing a multi-valued hyperparameter with the same number of values", {
  bk <- batch_kernel(se_kernel(length_scale = c(1, 2, 3)), batch_size = 3)
  bk2 <- kupdate(bk, length_scale = c(9, 8, 7))
  expect_equal(unname(sort(get_params(bk2)[grepl("length_scale", names(get_params(bk2)))])), c(7, 8, 9), tolerance = 1e-8)
})

test_that("kupdate() still allows a wrong length > 1 through (caught later, at evaluate())", {
  bk <- batch_kernel(se_kernel(length_scale = c(1, 2, 3)), batch_size = 3)
  bk2 <- kupdate(bk, length_scale = c(1, 2)) # no error here -- evaluate() re-validates the axis contract
  expect_error(evaluate(bk2, c(0, 1)), "batch_size")
})
