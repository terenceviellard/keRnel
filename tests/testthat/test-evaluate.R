test_that("evaluate() defaults x2 to x1 (auto-covariance)", {
  k <- se_kernel(length_scale = 1)
  expect_equal(evaluate(k, c(0, 1, 2)), evaluate(k, c(0, 1, 2), c(0, 1, 2)))
})

test_that("evaluate() accepts bare vectors and matrices interchangeably", {
  k <- se_kernel(length_scale = 1)
  x <- c(0, 1, 2)
  expect_equal(evaluate(k, x), evaluate(k, matrix(x, ncol = 1)))
})

test_that("evaluate() returns the right shape for unequal x1/x2 lengths", {
  k <- se_kernel(length_scale = 1)
  result <- evaluate(k, c(0, 1, 2), c(0, 1))
  expect_equal(dim(result), c(3, 2))
})

test_that("evaluate() propagates NA to the corresponding row/column", {
  k <- se_kernel(length_scale = 1)
  result <- evaluate(k, c(0, NA, 2), c(0, 1))
  expect_true(all(is.na(result[2, ])))
  expect_false(anyNA(result[1, ]))
  expect_false(anyNA(result[3, ]))
})

test_that("evaluate() propagates NA symmetrically for x2", {
  k <- se_kernel(length_scale = 1)
  result <- evaluate(k, c(0, 1), c(0, NA, 2))
  expect_true(all(is.na(result[, 2])))
  expect_false(anyNA(result[, 1]))
  expect_false(anyNA(result[, 3]))
})

test_that("evaluate() treats Inf like NA rather than producing a silent NaN", {
  k <- se_kernel(length_scale = 1)
  result <- evaluate(k, c(0, Inf, 2), c(0, 1))
  expect_true(all(is.na(result[2, ])))
  expect_false(anyNA(result[1, ]))
  expect_false(anyNA(result[3, ]))
  expect_false(any(is.nan(result)))
})

test_that("evaluate() treats -Inf like NA too, and respects na_action = \"mask\"", {
  k <- se_kernel(length_scale = 1)
  result <- evaluate(k, c(0, -Inf, 2), na_action = "mask")
  expect_false(anyNA(result))
  expect_false(any(is.nan(result)))
})

test_that("evaluate() dispatches via UseMethod on the kernel's class", {
  expect_true(is.function(getS3method("evaluate", "default")))
})

# na_action = "mask": keep the covariance matrix invertible despite missing
# data, by treating a missing point as an independent observation of
# variance `mask_variance` -- see dev/ameliorations/design/01-architecture.Rmd, "na_action".

test_that("na_action defaults to 'propagate' (unchanged behaviour)", {
  k <- se_kernel(length_scale = 1)
  result <- evaluate(k, c(0, NA, 2))
  expect_true(anyNA(result))
})

test_that("na_action = 'mask' produces no NA in the result", {
  k <- se_kernel(length_scale = 1)
  result <- evaluate(k, c(0, NA, 2), na_action = "mask")
  expect_false(anyNA(result))
})

test_that("na_action = 'mask' zeroes off-diagonal entries of a masked point", {
  k <- se_kernel(length_scale = 1)
  result <- evaluate(k, c(0, NA, 2), na_action = "mask")
  expect_equal(result[2, -2], c(0, 0))
  expect_equal(result[-2, 2], c(0, 0))
})

test_that("na_action = 'mask' sets the diagonal of a masked point to mask_variance", {
  k <- se_kernel(length_scale = 1)
  result <- evaluate(k, c(0, NA, 2), na_action = "mask")
  expect_equal(result[2, 2], 1)

  result2 <- evaluate(k, c(0, NA, 2), na_action = "mask", mask_variance = 3.5)
  expect_equal(result2[2, 2], 3.5)
})

test_that("na_action = 'mask' leaves non-missing entries unaffected", {
  k <- se_kernel(length_scale = 1)
  result <- evaluate(k, c(0, NA, 2), na_action = "mask")
  expect_equal(result[1, 1], 1)
  expect_equal(result[c(1, 3), c(1, 3)], evaluate(k, c(0, 2)))
})

test_that("na_action = 'mask' does not preserve a diagonal for cross-covariance", {
  k <- se_kernel(length_scale = 1)
  # x1 != x2: no diagonal to preserve, masked rows/cols are simply zeroed.
  result <- evaluate(k, c(0, NA, 2), c(0, 1), na_action = "mask")
  expect_equal(result[2, ], c(0, 0))
})

test_that("na_action = 'mask' with diagonal_engine sets mask_variance at masked indices", {
  k <- se_kernel(length_scale = 1, engine = diagonal_engine())
  result <- evaluate(k, c(0, NA, 2), c(0, 1, 2), na_action = "mask", mask_variance = 2)
  expect_equal(result[2], 2)
  expect_false(anyNA(result))
})

test_that("na_action rejects an invalid value", {
  k <- se_kernel(length_scale = 1)
  expect_error(evaluate(k, c(0, 1), na_action = "drop"))
})

test_that("na_action = 'mask' is forwarded through composite kernels", {
  k <- se_kernel(length_scale = 1) + white_noise_kernel(noise = 0.1)
  result <- evaluate(k, c(0, NA, 2), na_action = "mask")
  expect_false(anyNA(result))
})

test_that("na_action = 'mask' is forwarded through wrapper kernels", {
  k <- exp_kernel(se_kernel(length_scale = 1))
  result <- evaluate(k, c(0, NA, 2), na_action = "mask")
  expect_false(anyNA(result))
})

# Memory warning guard -- see dev/ameliorations/design/01-architecture.Rmd, "Garde-fou
# memoire".

test_that("evaluate() warns when the result would exceed memory_warning_threshold", {
  k <- se_kernel(length_scale = 1)
  expect_warning(evaluate(k, c(0, 1, 2), memory_warning_threshold = 2), "dense")
})

test_that("evaluate() does not warn below memory_warning_threshold", {
  k <- se_kernel(length_scale = 1)
  expect_no_warning(evaluate(k, c(0, 1, 2), memory_warning_threshold = 100))
})

test_that("memory_warning_threshold = Inf disables the warning entirely", {
  k <- se_kernel(length_scale = 1)
  expect_no_warning(evaluate(k, c(0, 1, 2), memory_warning_threshold = Inf))
})

test_that("evaluate() does not warn at the default threshold for small inputs", {
  k <- se_kernel(length_scale = 1)
  expect_no_warning(evaluate(k, c(0, 1, 2)))
})

# Multi-valued hyperparameter evaluated directly (finding A09/C06): only
# meaningful under a batching/blocking wrapper, which slices it down before
# evaluate.default() ever sees it. Evaluated directly, R would silently
# recycle the vector across the covariance matrix.

test_that("evaluate() rejects a base kernel with a multi-valued hyperparameter evaluated directly", {
  k <- se_kernel(length_scale = c(1, 3))
  expect_error(evaluate(k, c(0, 1, 2)), "more than 1 value")
})

test_that("evaluate() rejects a multi-valued hyperparameter set via kupdate() too", {
  k <- kupdate(se_kernel(length_scale = 1), length_scale = c(1, 3))
  expect_error(evaluate(k, c(0, 1, 2)), "more than 1 value")
})

test_that("evaluate() allows a block_compatible kernel a pair of values directly", {
  k <- feature_kernel(length_scale = c(1, 2), length_scale_u = 1, variance = c(1, 1))
  expect_no_error(evaluate(k, c(0, 1, 2)))
})

test_that("evaluate() still rejects more than a pair even for a block_compatible kernel", {
  k <- feature_kernel(length_scale = c(1, 2, 3), length_scale_u = 1, variance = 1)
  expect_error(evaluate(k, c(0, 1, 2)), "more than 2 value")
})

test_that("evaluate() does not flag a scalar hyperparameter", {
  k <- se_kernel(length_scale = 1)
  expect_no_error(evaluate(k, c(0, 1, 2)))
})

test_that("evaluate() through batch_kernel()/block_kernel()/input_specific_param_kernel() is unaffected", {
  expect_no_error(evaluate(batch_kernel(se_kernel(length_scale = c(1, 2, 3)), batch_size = 3), c(0, 1)))
  expect_no_error(evaluate(
    block_kernel(feature_kernel(length_scale = c(1, 2), length_scale_u = 1, variance = c(1, 1)), nb_blocks = 2),
    c(0, 1)
  ))
  expect_no_error(evaluate(
    input_specific_param_kernel(white_noise_kernel(noise = c(0.1, 0.5, 0.2)), input_size = 3),
    c(0, 1, 2)
  ))
})

test_that("evaluate_named() propagates rownames/names as result dimnames", {
  k <- se_kernel(length_scale = 1)
  x <- c(a = 0, b = 1, c = 2)
  result <- evaluate_named(k, x)
  expect_equal(unname(result), evaluate(k, x))
  expect_equal(rownames(result), c("a", "b", "c"))
  expect_equal(colnames(result), c("a", "b", "c"))
})

test_that("evaluate_named() uses x1's names for x2 when x2 is NULL, x2's own names otherwise", {
  k <- se_kernel(length_scale = 1)
  x1 <- matrix(c(0, 1), ncol = 1, dimnames = list(c("p1", "p2"), NULL))
  x2 <- matrix(c(2, 3, 4), ncol = 1, dimnames = list(c("q1", "q2", "q3"), NULL))
  result <- evaluate_named(k, x1, x2)
  expect_equal(rownames(result), c("p1", "p2"))
  expect_equal(colnames(result), c("q1", "q2", "q3"))
})

test_that("evaluate_named() returns the result unchanged when neither side has names", {
  k <- se_kernel(length_scale = 1)
  x <- c(0, 1, 2)
  expect_identical(evaluate_named(k, x), evaluate(k, x))
})

test_that("evaluate_named() names a diagonal_engine vector result", {
  k <- se_kernel(length_scale = 1, engine = diagonal_engine())
  x <- c(a = 0, b = 1, c = 2)
  result <- evaluate_named(k, x, x)
  expect_equal(names(result), c("a", "b", "c"))
})
