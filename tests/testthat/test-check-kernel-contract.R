# check_param_names(): diagnostic, never errors/warns on its own (finding #6,
# Option C -- a shared name is often intentional).

test_that("check_param_names() reports a name shared across sub-kernels", {
  k <- se_kernel(length_scale = 1) + periodic_kernel(length_scale = 2, period = 3)
  expect_equal(check_param_names(k), "length_scale")
})

test_that("check_param_names() returns character(0) when there is no collision", {
  k <- se_kernel(length_scale = 1) + linear_kernel(slope_var = 1)
  expect_equal(check_param_names(k), character(0))
})

test_that("check_param_names() does not error or warn", {
  k <- se_kernel(length_scale = 1) + se_kernel(length_scale = 2)
  expect_no_error(check_param_names(k))
  expect_no_warning(check_param_names(k))
})

# check_kernel_contract(): opt-in runtime diagnostics, complementary to
# validate_kernel().

test_that("check_kernel_contract() passes every check on a well-behaved kernel", {
  res <- check_kernel_contract(se_kernel(length_scale = 1))
  expect_true(res$symmetric)
  expect_true(res$positive_semidefinite)
  expect_true(res$self_covariance_consistent)
  expect_true(res$kupdate_roundtrip)
  expect_equal(res$shared_param_names, character(0))
})

test_that("check_kernel_contract() passes on a composite kernel with a legitimately shared name", {
  k <- se_kernel(length_scale = 1) + periodic_kernel(length_scale = 1, period = 3)
  res <- check_kernel_contract(k)
  expect_true(res$kupdate_roundtrip) # both occurrences share the same value (1): safe round-trip
  expect_equal(res$shared_param_names, "length_scale")
})

test_that("check_kernel_contract() reports NA for kupdate_roundtrip when a shared name has conflicting values", {
  k <- se_kernel(length_scale = 1) + se_kernel(length_scale = 2)
  res <- check_kernel_contract(k)
  expect_true(is.na(res$kupdate_roundtrip))
})

test_that("check_kernel_contract() reports NA for matrix-only checks under diagonal_engine()", {
  k <- se_kernel(length_scale = 1, engine = diagonal_engine())
  res <- check_kernel_contract(k, X = c(0, 1, 2))
  expect_true(is.na(res$symmetric))
  expect_true(is.na(res$positive_semidefinite))
})

test_that("check_kernel_contract() rejects a non-kernel argument", {
  expect_error(check_kernel_contract(1), "kernel")
})

test_that("print.kernel_contract_check() prints a pass/fail summary", {
  res <- check_kernel_contract(se_kernel(length_scale = 1))
  expect_output(print(res), "check_kernel_contract")
  expect_output(print(res), "symmetric")
})
