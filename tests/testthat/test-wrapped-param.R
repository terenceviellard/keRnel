test_that("wrapped_param round-trips the constrained value", {
  p <- wrapped_param(2.5)
  expect_equal(unwrap_param(p), 2.5, tolerance = 1e-8)
})

test_that("wrapped_param defaults to log_exp_parametrisation", {
  p <- wrapped_param(1)
  expect_s3_class(p$parametrisation, "log_exp_parametrisation")
})

test_that("wrapped_param respects a custom parametrisation", {
  p <- wrapped_param(-2, parametrisation = identity_parametrisation())
  expect_equal(unwrap_param(p), -2)
})

test_that("wrapped_param defaults to not frozen, and frozen can be set explicitly", {
  expect_false(wrapped_param(1)$frozen)
  expect_true(wrapped_param(1, frozen = TRUE)$frozen)
})

test_that("invalid values are rejected at construction (validation lives in wrap())", {
  expect_error(wrapped_param(-1), "strictly positive")
})

test_that("wrapped_param stores the value in the unconstrained space", {
  p <- wrapped_param(1) # log_exp_parametrisation: wrap(1) == log(1) == 0
  expect_equal(p$wrapped, 0)
})

test_that("format.wrapped_param shows the constrained value and frozen status", {
  expect_match(format(wrapped_param(2)), "^2")
  expect_no_match(format(wrapped_param(2)), "frozen")
  expect_match(format(wrapped_param(2, frozen = TRUE)), "frozen")
})

test_that("print.wrapped_param prints the formatted value and returns its argument invisibly", {
  p <- wrapped_param(2)
  expect_output(print(p), "^2")
  expect_identical(withVisible(print(p))$visible, FALSE)
})

test_that("format.wrapped_param collapses a vector-valued parameter to a single string", {
  # Regression test: a naive format(value) on a length > 1 value returns one
  # string per element, which broke format.kernel()'s vapply() once
  # vector-valued hyperparameters were introduced (ard_kernel(), batch_kernel()).
  p <- wrapped_param(c(1, 2, 3))
  formatted <- format(p)
  expect_length(formatted, 1)
  expect_equal(formatted, "[1, 2, 3]")
})

test_that("format.wrapped_param shows frozen status for a vector-valued parameter", {
  p <- wrapped_param(c(1, 2), frozen = TRUE)
  expect_match(format(p), "frozen")
  expect_length(format(p), 1)
})
