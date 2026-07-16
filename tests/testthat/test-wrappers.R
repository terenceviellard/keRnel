# exp_kernel: k(x1, x2) = exp(k_inner(x1, x2))

test_that("exp_kernel applies exp() to its inner kernel's evaluate() result", {
  inner <- se_kernel(length_scale = 1)
  k <- exp_kernel(inner)
  expect_equal(evaluate(k, c(0, 1, 2)), exp(evaluate(inner, c(0, 1, 2))))
})

test_that("exp_kernel's class vector reflects its category hierarchy", {
  k <- exp_kernel(se_kernel(length_scale = 1))
  expect_equal(class(k), c("exp_kernel", "wrapper_kernel", "kernel"))
})

test_that("exp_kernel coerces a bare number to a constant_kernel", {
  k <- exp_kernel(2)
  expect_true(inherits(k$subkernels[[1]], "constant_kernel"))
  expect_equal(evaluate(k, c(0, 1))[1, 1], exp(2))
})

test_that("exp_kernel's format() wraps the inner kernel's format()", {
  k <- exp_kernel(se_kernel(length_scale = 1))
  expect_match(format(k), "^Exp\\(.*\\)$")
})

# log_kernel: k(x1, x2) = log(k_inner(x1, x2))

test_that("log_kernel applies log() to its inner kernel's evaluate() result", {
  inner <- variance_kernel(variance = 3)
  k <- log_kernel(inner)
  expect_equal(evaluate(k, c(0, 1))[1, 1], log(3))
})

test_that("exp_kernel and log_kernel are inverses of each other", {
  inner <- se_kernel(length_scale = 1)
  k <- log_kernel(exp_kernel(inner))
  expect_equal(evaluate(k, c(0, 1, 2)), evaluate(inner, c(0, 1, 2)), tolerance = 1e-10)
})

test_that("log_kernel's format() wraps the inner kernel's format()", {
  k <- log_kernel(variance_kernel(variance = 1))
  expect_match(format(k), "^Log\\(.*\\)$")
})

# get_params()/kupdate()/freeze() recurse through wrapper kernels for free

test_that("get_params() collects hyperparameters through a wrapper kernel", {
  k <- exp_kernel(se_kernel(length_scale = 4))
  expect_equal(unname(get_params(k)), 4, tolerance = 1e-10)
})

test_that("kupdate() updates a hyperparameter inside a wrapper kernel", {
  k <- log_kernel(se_kernel(length_scale = 1))
  k2 <- kupdate(k, length_scale = 7)
  expect_equal(unname(get_params(k2)), 7, tolerance = 1e-10)
})

# Validation

test_that("validate_kernel() rejects a hand-built wrapper_kernel with a non-kernel subkernel", {
  broken <- structure(
    list(subkernels = list(4)),
    class = c("broken_kernel", "wrapper_kernel", "kernel")
  )
  expect_error(validate_kernel(broken), "must itself be a `kernel`")
})

test_that("validate_kernel() rejects a hand-built wrapper_kernel with the wrong number of subkernels", {
  broken <- structure(
    list(subkernels = list(se_kernel(length_scale = 1), se_kernel(length_scale = 1))),
    class = c("broken_kernel", "wrapper_kernel", "kernel")
  )
  expect_error(validate_kernel(broken), "exactly 1")
})

# active_dims_kernel: restricts the inner kernel to a subset of columns,
# selected before delegating to evaluate().

test_that("active_dims_kernel only exposes the selected column(s) to the inner kernel", {
  k <- active_dims_kernel(se_kernel(length_scale = 1), active_dims = 2)
  X <- matrix(c(0, 10, 1, 11, 2, 12), ncol = 2, byrow = TRUE)
  expect_equal(evaluate(k, X), evaluate(se_kernel(length_scale = 1), c(10, 11, 12)))
})

test_that("active_dims_kernel supports more than one selected dimension", {
  k <- active_dims_kernel(se_kernel(length_scale = 1), active_dims = c(1, 3))
  X <- matrix(c(0, 99, 1, 1, 99, 2, 2, 99, 3), ncol = 3, byrow = TRUE)
  expect_equal(evaluate(k, X), evaluate(se_kernel(length_scale = 1), X[, c(1, 3)]))
})

test_that("active_dims_kernel handles distinct x1/x2", {
  k <- active_dims_kernel(linear_kernel(slope_var = 1), active_dims = 1)
  X1 <- matrix(c(1, 9, 2, 9), ncol = 2, byrow = TRUE)
  X2 <- matrix(c(3, 9, 4, 9), ncol = 2, byrow = TRUE)
  expect_equal(evaluate(k, X1, X2), evaluate(linear_kernel(slope_var = 1), X1[, 1], X2[, 1]))
})

test_that("active_dims_kernel errors when active_dims is out of bounds", {
  k <- active_dims_kernel(se_kernel(length_scale = 1), active_dims = 3)
  X <- matrix(c(0, 0, 1, 1), ncol = 2, byrow = TRUE)
  expect_error(evaluate(k, X), "out of bounds")
})

test_that("active_dims_kernel rejects non-positive active_dims at construction", {
  expect_error(active_dims_kernel(se_kernel(length_scale = 1), active_dims = 0), "positive integers")
  expect_error(active_dims_kernel(se_kernel(length_scale = 1), active_dims = -1), "positive integers")
})

test_that("active_dims_kernel's format() shows the selected dimensions", {
  k <- active_dims_kernel(se_kernel(length_scale = 1), active_dims = c(2, 3))
  expect_match(format(k), "\\[dims=2,3\\]")
})

test_that("get_params() and kupdate() pass through active_dims_kernel", {
  k <- active_dims_kernel(se_kernel(length_scale = 2), active_dims = 1)
  expect_equal(unname(get_params(k)), 2, tolerance = 1e-10)
  k2 <- kupdate(k, length_scale = 9)
  expect_equal(unname(get_params(k2)), 9, tolerance = 1e-10)
})

# ard_kernel: gives every input dimension its own length_scale, by
# rescaling the inputs before delegating to a stationary kernel built (and
# frozen) internally -- see R/wrappers.R for why the inner kernel's own
# length_scale can't accidentally be left trainable and misused.

test_that("ard_kernel matches manually rescaling inputs by length_scales", {
  k <- ard_kernel(se_kernel, length_scales = c(2, 3))
  X <- matrix(c(0, 0, 2, 3, -2, 6), ncol = 2, byrow = TRUE)
  rescaled <- sweep(X, 2, c(2, 3), "/")
  expect_equal(evaluate(k, X), evaluate(se_kernel(length_scale = 1), rescaled))
})

test_that("ard_kernel forwards extra constructor arguments via ...", {
  k <- ard_kernel(periodic_kernel, length_scales = c(1, 2), period = 5)
  expect_equal(unwrap_param(k$subkernels[[1]]$period), 5)
})

test_that("ard_kernel's inner length_scale is fixed at 1 and frozen", {
  k <- ard_kernel(se_kernel, length_scales = c(1, 2, 3))
  inner <- k$subkernels[[1]]
  expect_equal(unwrap_param(inner$length_scale), 1)
  expect_true(inner$length_scale$frozen)
})

test_that("ard_kernel's length_scales are the only trainable parameters by default", {
  k <- ard_kernel(se_kernel, length_scales = c(1, 2, 3))
  expect_equal(unname(get_trainable_params(k)), c(1, 2, 3), tolerance = 1e-10)
  expect_length(get_params(k), 4) # length_scales (3) + frozen inner length_scale (1)
})

test_that("ard_kernel rejects non-positive length_scales", {
  expect_error(ard_kernel(se_kernel, length_scales = c(1, -1)), "strictly positive")
})

test_that("ard_kernel errors when input dimensionality does not match length_scales", {
  k <- ard_kernel(se_kernel, length_scales = c(1, 2, 3))
  X <- matrix(rnorm(4), ncol = 2)
  expect_error(evaluate(k, X), "length_scales")
})

test_that("kupdate() updates ard_kernel's length_scales", {
  k <- ard_kernel(se_kernel, length_scales = c(1, 2))
  k2 <- kupdate(k, length_scales = c(4, 5))
  X <- matrix(c(0, 0, 1, 1), ncol = 2, byrow = TRUE)
  expect_equal(evaluate(k2, X), evaluate(se_kernel(length_scale = 1), sweep(X, 2, c(4, 5), "/")))
})

test_that("ard_kernel's format() shows the inner kernel and length_scales", {
  k <- ard_kernel(se_kernel, length_scales = c(1, 2))
  expect_match(format(k), "^ARD\\(")
})

# slice_kernel(): the "Share" primitive backing batch_kernel() and
# input_specific_param_kernel() -- only the hyperparameters named in `axes`
# (the contract resolved by .resolve_axes()) are sliced at index b; a
# length-1 (shared) value, or a hyperparameter not named in `axes`, is
# left untouched.

test_that("slice_kernel() slices an axis hyperparameter at the given index", {
  k <- se_kernel(length_scale = c(1, 2, 3))
  k2 <- slice_kernel(k, 2, axes = "length_scale")
  expect_equal(unwrap_param(k2$length_scale), 2)
})

test_that("slice_kernel() leaves a shared (length-1) hyperparameter untouched", {
  k <- se_kernel(length_scale = 5)
  k2 <- slice_kernel(k, 3, axes = "length_scale")
  expect_equal(unwrap_param(k2$length_scale), 5)
})

test_that("slice_kernel() leaves a hyperparameter not named in axes untouched", {
  k <- se_kernel(length_scale = c(1, 2)) + white_noise_kernel(noise = c(10, 20))
  k2 <- slice_kernel(k, 2, axes = "noise")
  expect_equal(unwrap_param(k2$subkernels[[1]]$length_scale), c(1, 2))
  expect_equal(unwrap_param(k2$subkernels[[2]]$noise), 20)
})

test_that("slice_kernel() recurses through $subkernels", {
  k <- se_kernel(length_scale = c(1, 2)) + white_noise_kernel(noise = c(10, 20))
  k2 <- slice_kernel(k, 2, axes = c("length_scale", "noise"))
  expect_equal(unwrap_param(k2$subkernels[[1]]$length_scale), 2)
  expect_equal(unwrap_param(k2$subkernels[[2]]$noise), 20)
})

# batch_kernel: each batch element may have its own hyperparameter values.
# Which hyperparameters vary per batch element is an explicit contract
# (.resolve_axes()): auto-detected from value lengths by default (and then
# strictly validated against batch_size), or named via `batch_axes`.

test_that("batch_kernel produces one covariance matrix per batch element", {
  k <- batch_kernel(se_kernel(length_scale = c(1, 2, 3)), batch_size = 3)
  result <- evaluate(k, c(0, 1, 2))
  expect_equal(dim(result), c(3, 3, 3))
  for (b in seq_len(3)) {
    expect_equal(result[, , b], evaluate(se_kernel(length_scale = b), c(0, 1, 2)))
  }
})

test_that("batch_kernel shares a scalar hyperparameter across every batch element", {
  k <- batch_kernel(se_kernel(length_scale = 2), batch_size = 4)
  result <- evaluate(k, c(0, 1))
  reference <- evaluate(se_kernel(length_scale = 2), c(0, 1))
  for (b in seq_len(4)) {
    expect_equal(result[, , b], reference)
  }
})

test_that("batch_kernel with batch_over_inputs uses a different x1 per batch element", {
  k <- batch_kernel(se_kernel(length_scale = 1), batch_size = 2, batch_over_inputs = TRUE)
  x1 <- array(0, dim = c(2, 3, 1))
  x1[1, , ] <- c(0, 1, 2)
  x1[2, , ] <- c(10, 11, 12)
  result <- evaluate(k, x1)
  expect_equal(result[, , 1], evaluate(se_kernel(length_scale = 1), c(0, 1, 2)))
  expect_equal(result[, , 2], evaluate(se_kernel(length_scale = 1), c(10, 11, 12)))
})

test_that("batch_kernel rejects a non-3D input when batch_over_inputs = TRUE", {
  k <- batch_kernel(se_kernel(length_scale = 1), batch_size = 2, batch_over_inputs = TRUE)
  expect_error(evaluate(k, c(0, 1, 2)), "3-dimensional")
})

test_that("batch_kernel rejects a non-positive batch_size", {
  expect_error(batch_kernel(se_kernel(length_scale = 1), batch_size = 0), "positive integer")
})

test_that("batch_kernel preserves the (1, 1, batch_size) array shape for a single input point (dense engine)", {
  k <- batch_kernel(se_kernel(length_scale = c(1, 2, 3)), batch_size = 3)
  result <- evaluate(k, matrix(5, nrow = 1))
  expect_equal(dim(result), c(1, 1, 3))
  for (b in 1:3) {
    expect_equal(as.numeric(result[, , b]), as.numeric(evaluate(se_kernel(length_scale = b), matrix(5, nrow = 1))))
  }
})

test_that("batch_kernel preserves the (1, batch_size) matrix shape for a single input point (diagonal engine)", {
  k <- batch_kernel(se_kernel(length_scale = c(1, 2), engine = diagonal_engine()), batch_size = 2)
  result <- evaluate(k, matrix(5, nrow = 1), matrix(5, nrow = 1))
  expect_equal(dim(result), c(1, 2))
})

test_that("block_diag_kernel() with nb_blocks = 1 no longer crashes on a single input point", {
  k <- block_diag_kernel(se_kernel(length_scale = 2), nb_blocks = 1)
  expect_equal(evaluate(k, matrix(5, nrow = 1)), evaluate(se_kernel(length_scale = 2), matrix(5, nrow = 1)))
})

# The plan-based fast path (R/plan.R, precompute()/eval_plan()) routes
# evaluate.batch_kernel() through a cached distance plan whenever
# batch_over_inputs = FALSE -- these tests check it reproduces exactly what
# evaluating each slice directly (the pre-existing, always-correct path)
# would produce, for every case the plan must fall back on or combine.

test_that("batch_kernel's plan-based path matches direct evaluation for a constant_kernel (flat node) sibling", {
  inner <- se_kernel(length_scale = c(1, 2)) + constant_kernel(value = 3)
  k <- batch_kernel(inner, batch_size = 2)
  result <- evaluate(k, c(0, 1, 2))
  for (b in 1:2) {
    reference <- evaluate(se_kernel(length_scale = b) + constant_kernel(value = 3), c(0, 1, 2))
    expect_equal(result[, , b], reference)
  }
})

test_that("batch_kernel's plan-based path falls back correctly for an affine_kernel sibling (HP-dependent distance)", {
  inner <- se_kernel(length_scale = c(1, 2, 3)) + affine_kernel(slope_var = 2, offset = 0.5)
  k <- batch_kernel(inner, batch_size = 3)
  result <- evaluate(k, c(0, 1, 2))
  for (b in 1:3) {
    reference <- evaluate(se_kernel(length_scale = b) + affine_kernel(slope_var = 2, offset = 0.5), c(0, 1, 2))
    expect_equal(result[, , b], reference)
  }
})

test_that("batch_kernel's plan-based path handles a diagonal_engine leaf", {
  inner <- se_kernel(length_scale = c(1, 2), engine = diagonal_engine())
  k <- batch_kernel(inner, batch_size = 2)
  result <- evaluate(k, c(0, 1, 2), c(0, 1, 2))
  expect_equal(dim(result), c(3, 2))
  for (b in 1:2) {
    reference <- evaluate(se_kernel(length_scale = b, engine = diagonal_engine()), c(0, 1, 2), c(0, 1, 2))
    expect_equal(result[, b], reference)
  }
})

test_that("batch_kernel's plan-based path masks NA per-leaf, matching direct evaluation of a sum of two leaves", {
  inner <- se_kernel(length_scale = c(1, 2)) + white_noise_kernel(noise = 0.1)
  k <- batch_kernel(inner, batch_size = 2)
  x <- c(0, NA, 2)
  for (na_mode in c("propagate", "mask")) {
    result <- evaluate(k, x, na_action = na_mode)
    for (b in 1:2) {
      reference <- evaluate(se_kernel(length_scale = b) + white_noise_kernel(noise = 0.1), x, na_action = na_mode)
      expect_equal(result[, , b], reference)
    }
  }
})

test_that("batch_kernel's plan-based path reuses one plan correctly across multi-axis combos", {
  k3 <- variance_kernel(variance = 1:6) * se_kernel(length_scale = 1)
  mk <- batch_kernel(k3, batch_size = c(K = 3, O = 2), batch_axes = list(variance = c("K", "O")))
  result <- evaluate(mk, c(0, 1, 2))
  expect_equal(dim(result), c(3, 3, 3, 2))
  combo_value <- matrix(1:6, nrow = 3, ncol = 2) # column-major: K fastest-varying
  for (kk in 1:3) {
    for (oo in 1:2) {
      reference <- evaluate(variance_kernel(variance = combo_value[kk, oo]) * se_kernel(length_scale = 1), c(0, 1, 2))
      expect_equal(result[, , kk, oo], reference)
    }
  }
})

test_that("get_params() and kupdate() pass through a batch_kernel", {
  k <- batch_kernel(se_kernel(length_scale = c(1, 2)), batch_size = 2)
  expect_equal(unname(get_params(k)), c(1, 2), tolerance = 1e-10)
  k2 <- kupdate(k, length_scale = c(9, 8))
  expect_equal(unname(get_params(k2)), c(9, 8), tolerance = 1e-10)
})

test_that("a nested batch_kernel (double batching) still works through $subkernels", {
  # Verifies that nested batching composes correctly without special-casing:
  # freeze()/get_params()/kupdate() never rely on a separate mask tree that
  # could go out of sync; they read `frozen`/wrapped values directly off each
  # wrapped_param, recursing through $subkernels uniformly regardless of
  # nesting depth.
  inner <- batch_kernel(se_kernel(length_scale = c(1, 2)), batch_size = 2)
  outer <- batch_kernel(inner, batch_size = 2)

  expect_equal(unname(get_params(outer)), c(1, 2), tolerance = 1e-10)

  outer_frozen <- freeze(outer, "length_scale")
  expect_length(get_trainable_params(outer_frozen), 0)

  outer2 <- kupdate(outer, length_scale = c(5, 6))
  expect_equal(unname(get_params(outer2)), c(5, 6), tolerance = 1e-10)
})

test_that("batch_kernel's format() shows just the inner kernel", {
  k <- batch_kernel(se_kernel(length_scale = c(1, 2)), batch_size = 2)
  expect_equal(format(k), format(k$subkernels[[1]]))
})

# input_specific_param_kernel: each point of x1 gets its own hyperparameter
# value, under the same explicit axis contract as batch_kernel
# (auto-detected by default, or named via `param_axes`).

test_that("input_specific_param_kernel gives each point of x1 its own hyperparameter value", {
  k <- input_specific_param_kernel(
    white_noise_kernel(noise = c(0.1, 0.5, 0.2)),
    input_size = 3
  )
  result <- evaluate(k, c(0, 1, 2))
  expect_equal(diag(result), c(0.1, 0.5, 0.2))
  expect_equal(result[upper.tri(result)], rep(0, 3))
})

test_that("input_specific_param_kernel shares a scalar hyperparameter across all points", {
  k <- input_specific_param_kernel(se_kernel(length_scale = 1.5), input_size = 3)
  result <- evaluate(k, c(0, 1, 2))
  expect_equal(result, evaluate(se_kernel(length_scale = 1.5), c(0, 1, 2)), tolerance = 1e-10)
})

test_that("input_specific_param_kernel errors when nrow(x1) does not match input_size", {
  k <- input_specific_param_kernel(se_kernel(length_scale = 1), input_size = 3)
  expect_error(evaluate(k, c(0, 1)), "input_size")
})

test_that("input_specific_param_kernel rejects a non-positive input_size", {
  expect_error(input_specific_param_kernel(se_kernel(length_scale = 1), input_size = 0), "positive integer")
})

test_that("get_params() and kupdate() pass through an input_specific_param_kernel", {
  k <- input_specific_param_kernel(white_noise_kernel(noise = c(1, 2, 3)), input_size = 3)
  expect_equal(unname(get_params(k)), c(1, 2, 3), tolerance = 1e-10)
  k2 <- kupdate(k, noise = c(4, 5, 6))
  expect_equal(unname(get_params(k2)), c(4, 5, 6), tolerance = 1e-10)
})

test_that("input_specific_param_kernel's format() shows just the inner kernel", {
  k <- input_specific_param_kernel(se_kernel(length_scale = 1), input_size = 2)
  expect_equal(format(k), format(k$subkernels[[1]]))
})

# block_diag_kernel: reuses batch_kernel(), then assembles block-diagonally
# (off-diagonal blocks are 0).

test_that("block_diag_kernel places each batch element's matrix on the diagonal", {
  k <- block_diag_kernel(se_kernel(length_scale = c(1, 2)), nb_blocks = 2)
  x <- c(0, 1)
  result <- evaluate(k, x)

  expect_equal(dim(result), c(4, 4))
  expect_equal(result[1:2, 1:2], evaluate(se_kernel(length_scale = 1), x))
  expect_equal(result[3:4, 3:4], evaluate(se_kernel(length_scale = 2), x))
})

test_that("block_diag_kernel's off-diagonal blocks are all zero", {
  k <- block_diag_kernel(se_kernel(length_scale = c(1, 2)), nb_blocks = 2)
  result <- evaluate(k, c(0, 1))
  expect_equal(result[1:2, 3:4], matrix(0, 2, 2))
  expect_equal(result[3:4, 1:2], matrix(0, 2, 2))
})

test_that("block_diag_kernel's format() shows BlockDiag + inner", {
  k <- block_diag_kernel(se_kernel(length_scale = c(1, 2)), nb_blocks = 2)
  expect_match(format(k), "^BlockDiag")
})

test_that("get_params() and kupdate() pass through a block_diag_kernel", {
  k <- block_diag_kernel(se_kernel(length_scale = c(1, 2)), nb_blocks = 2)
  expect_equal(unname(get_params(k)), c(1, 2), tolerance = 1e-10)
  k2 <- kupdate(k, length_scale = c(7, 8))
  result <- evaluate(k2, c(0, 1))
  expect_equal(result[3:4, 3:4], evaluate(se_kernel(length_scale = 8), c(0, 1)))
})

# pair_slice_kernel(): the block_kernel() counterpart of slice_kernel(),
# extracting a *pair* of elements rather than a single one.

test_that("pair_slice_kernel() extracts a pair of elements from an axis hyperparameter", {
  k <- se_kernel(length_scale = c(10, 20, 30))
  k2 <- pair_slice_kernel(k, 1, 3, axes = "length_scale")
  expect_equal(unwrap_param(k2$length_scale), c(10, 30))
})

test_that("pair_slice_kernel() leaves a shared (length-1) hyperparameter untouched", {
  k <- se_kernel(length_scale = 5)
  k2 <- pair_slice_kernel(k, 1, 2, axes = "length_scale")
  expect_equal(unwrap_param(k2$length_scale), 5)
})

# block_kernel: full block matrix, including cross-block terms, restricted
# to inner kernels explicitly marked block_compatible when hyperparameters
# vary per block (feature_kernel is the prototype -- see R/stationary-kernels.R).

test_that("block_kernel rejects a non-block_compatible inner with per-block hyperparameters", {
  expect_error(
    block_kernel(se_kernel(length_scale = c(1, 2)), nb_blocks = 2),
    "block_compatible"
  )
})

test_that("block_kernel accepts a non-block_compatible inner with only shared hyperparameters", {
  k <- block_kernel(se_kernel(length_scale = 1), nb_blocks = 2)
  result <- evaluate(k, c(0, 1))
  expect_equal(dim(result), c(4, 4))
  # Every block is identical (the whole point is tiling, no pairing needed).
  expect_equal(result[1:2, 1:2], result[3:4, 3:4])
  expect_equal(result[1:2, 3:4], result[1:2, 1:2])
})

test_that("block_kernel computes non-zero cross-block terms for a block_compatible inner", {
  inner <- feature_kernel(length_scale = c(1, 3), length_scale_u = 1, variance = c(1, 1))
  k <- block_kernel(inner, nb_blocks = 2)
  result <- evaluate(k, c(0, 1))
  expect_equal(dim(result), c(4, 4))
  expect_false(all(result[1:2, 3:4] == 0))
})

test_that("block_kernel agrees with a naive (non-symmetry-shortcut) reference", {
  inner <- feature_kernel(length_scale = c(1, 3), length_scale_u = 1, variance = c(2, 1))
  k <- block_kernel(inner, nb_blocks = 2)
  result <- evaluate(k, c(0, 1, 2))

  naive_blocks <- lapply(1:2, function(i) {
    lapply(1:2, function(j) {
      evaluate(pair_slice_kernel(inner, i, j, axes = c("length_scale", "variance")), c(0, 1, 2))
    })
  })
  naive <- do.call(rbind, lapply(naive_blocks, function(row) do.call(cbind, row)))

  expect_equal(result, naive, tolerance = 1e-10)
})

test_that("block_kernel's cross-block term is the transpose of its mirror block", {
  inner <- feature_kernel(length_scale = c(1, 3), length_scale_u = 1, variance = c(2, 1))
  k <- block_kernel(inner, nb_blocks = 2)
  result <- evaluate(k, c(0, 1))
  expect_equal(result[1:2, 3:4], t(result[3:4, 1:2]))
})

test_that("block_kernel rejects a non-positive nb_blocks", {
  expect_error(block_kernel(se_kernel(length_scale = 1), nb_blocks = 0), "positive integer")
})

test_that("block_kernel's format() shows Block + inner", {
  k <- block_kernel(se_kernel(length_scale = 1), nb_blocks = 2)
  expect_match(format(k), "^Block")
})

test_that("get_params() and kupdate() pass through a block_kernel", {
  inner <- feature_kernel(length_scale = c(1, 2), length_scale_u = 1, variance = c(1, 1))
  k <- block_kernel(inner, nb_blocks = 2)
  expect_equal(unname(get_params(k))[1:2], c(1, 2), tolerance = 1e-10)

  k2 <- kupdate(k, length_scale = c(5, 6))
  updated_inner <- kupdate(inner, length_scale = c(5, 6))
  expect_equal(
    evaluate(k2, c(0, 1)),
    evaluate(block_kernel(updated_inner, nb_blocks = 2), c(0, 1)),
    tolerance = 1e-10
  )
})

# The explicit axis contract (.resolve_axes()): which hyperparameters vary
# along a wrapper's batch/block/input axis is declared (or strictly
# auto-detected), validated at construction AND re-validated at evaluate()
# time -- never inferred per-slice from raw vector lengths. These are the
# regression tests for the silent-mis-slicing family of bugs (see
# dev/ameliorations/design/04-audit-bugs-fonctionnalites-wrappers.Rmd and
# dev/ameliorations/design/audit-multi-agents/).

test_that("batch_kernel errors at construction when a hyperparameter length does not match batch_size", {
  # Too few values: indexing past the end used to yield silent NAs.
  expect_error(
    batch_kernel(se_kernel(length_scale = c(1, 2, 3)), batch_size = 5),
    "batch_size"
  )
  # Too many values: the extra ones used to be silently ignored.
  expect_error(
    batch_kernel(se_kernel(length_scale = c(1, 2, 3)), batch_size = 2),
    "batch_size"
  )
})

test_that("batch_kernel errors when an explicitly named batch axis has the wrong length", {
  expect_error(
    batch_kernel(se_kernel(length_scale = c(1, 2, 3)), batch_size = 5,
                 batch_axes = "length_scale"),
    "length_scale"
  )
})

test_that("batch_kernel errors on a batch_axes name that does not exist", {
  expect_error(
    batch_kernel(se_kernel(length_scale = 1), batch_size = 2,
                 batch_axes = "lenght_scale"),
    "No such hyperparameter"
  )
})

test_that("evaluate() re-validates the batch axis contract after kupdate()", {
  bk <- batch_kernel(se_kernel(length_scale = c(1, 2, 3)), batch_size = 3)
  bk2 <- kupdate(bk, length_scale = c(1, 2)) # now length 2 != batch_size 3
  expect_error(evaluate(bk2, c(0, 1)), "batch_size")
})

test_that("ard_kernel's per-dimension length_scales are never mistaken for a batch axis", {
  # batch_size == length(length_scales) == 3: under pure length-based
  # inference this coincidence used to slice the per-*dimension* vector
  # into a meaningless per-batch scalar, silently losing ARD semantics.
  ard <- ard_kernel(se_kernel, length_scales = c(1, 2, 3))
  bk <- batch_kernel(ard, batch_size = 3) # auto-detect: nothing is batched
  X <- matrix(seq_len(12), ncol = 3)
  result <- evaluate(bk, X)
  expect_equal(dim(result), c(4, 4, 3))
  reference <- evaluate(ard, X)
  for (b in 1:3) {
    expect_equal(result[, , b], reference, tolerance = 1e-12)
  }
})

test_that("naming a non-sliceable hyperparameter as a batch axis errors", {
  ard <- ard_kernel(se_kernel, length_scales = c(1, 2, 3))
  expect_error(
    batch_kernel(ard, batch_size = 3, batch_axes = "length_scales"),
    "input dimension"
  )
})

test_that("nested batching with different sizes per level requires explicit axes", {
  # N.b. unlike kernax's equivalent test, noise values must be strictly
  # positive here: log_exp_parametrisation() rejects 0 at construction
  # (kernax silently wraps it to log(0) = -Inf).
  k <- se_kernel(length_scale = c(.5, 1.5)) + white_noise_kernel(noise = c(0.5, 1, 2, 3))
  # Auto-detection cannot decide which of the two multi-valued
  # hyperparameters belongs to which nesting level: loud error, pointing
  # at batch_axes.
  expect_error(batch_kernel(k, batch_size = 2), "batch_axes")
})

test_that("nested batch_kernel with explicit axes reproduces kernax's double-batch semantics", {
  # The R counterpart of kernax's test_double_batch_with_masks
  # (tests/test_kernel_wrappers.py): inner batch (size 2) varies
  # length_scale, outer batch (size 4) varies noise. Each level slices
  # only its own axis; result gains one trailing dimension per level.
  k <- se_kernel(length_scale = c(.5, 1.5)) + white_noise_kernel(noise = c(0.5, 1, 2, 3))
  bk <- batch_kernel(k, batch_size = 2, batch_axes = "length_scale")
  bbk <- batch_kernel(bk, batch_size = 4, batch_axes = "noise")

  result <- evaluate(bbk, c(1, 2, 3))
  expect_equal(dim(result), c(3, 3, 2, 4))

  ls_vals <- c(.5, 1.5)
  noise_vals <- c(0.5, 1, 2, 3)
  for (o in 1:4) {
    for (i in 1:2) {
      expected <- evaluate(
        se_kernel(length_scale = ls_vals[i]) + white_noise_kernel(noise = noise_vals[o]),
        c(1, 2, 3)
      )
      expect_equal(result[, , i, o], expected, tolerance = 1e-12)
    }
  }

  # Spot-checks against the same closed-form reference values as the
  # kernax test: diagonal = 1 + noise (independent of the inner batch),
  # off-diagonal = exp(-0.5 d^2 / l^2) (independent of the outer batch).
  expect_equal(result[1, 2, 1, 1], exp(-2), tolerance = 1e-12)   # l = .5, d = 1
  expect_equal(result[1, 3, 1, 1], exp(-8), tolerance = 1e-12)   # l = .5, d = 2
  expect_equal(result[1, 2, 2, 1], exp(-2 / 9), tolerance = 1e-12) # l = 1.5, d = 1
  expect_equal(diag(result[, , 1, 3]), rep(3, 3), tolerance = 1e-12) # 1 + noise, noise = 2
})

# Multi-axis batch_kernel: a single wrapper handling several *named* batch
# dimensions at once (`batch_size = c(name = size, ...)`), so that one
# hyperparameter can vary along the full Cartesian product of 2+ axes --
# something nested single-axis batch_kernel() calls cannot express, because
# each level collapses an axis hyperparameter to a scalar before the inner
# level sees it. Motivated by porting MIMOSA (dev/ameliorations/design/07): its
# shared_cluster_hps/shared_output_hps flags require e.g. a mean-kernel
# hyperparameter with one independent value per (cluster, output)
# combination, not per cluster or per output alone.

test_that("multi-axis batch_kernel produces one array dimension per named axis", {
  k <- variance_kernel(variance = 1:6) * se_kernel(length_scale = 1)
  mk <- batch_kernel(k, batch_size = c(K = 3, O = 2),
                      batch_axes = list(variance = c("K", "O")))
  result <- evaluate(mk, c(0, 1, 2))
  expect_equal(dim(result), c(3, 3, 3, 2))
})

test_that("multi-axis batch_kernel indexes a cross-axis hyperparameter in column-major order", {
  # axes = c("K", "O") with batch_size = c(K = 3, O = 2): K is the
  # fastest-varying axis, exactly like array(1:6, dim = c(3, 2)) -- flat
  # index 1 + (k - 1) + (o - 1) * 3.
  k <- variance_kernel(variance = 1:6) * se_kernel(length_scale = 1)
  mk <- batch_kernel(k, batch_size = c(K = 3, O = 2),
                      batch_axes = list(variance = c("K", "O")))
  result <- evaluate(mk, c(0, 1, 2))

  for (o in 1:2) {
    for (kk in 1:3) {
      flat <- kk + (o - 1) * 3
      expected <- evaluate(variance_kernel(variance = flat) * se_kernel(length_scale = 1), c(0, 1, 2))
      expect_equal(result[, , kk, o], expected, tolerance = 1e-12)
    }
  }
})

test_that("multi-axis batch_kernel broadcasts a hyperparameter declared on a single axis", {
  # `variance` varies only along K; every O slice for a given K must be
  # identical, and equal to the single-axis reference.
  k <- variance_kernel(variance = c(10, 20, 30)) * se_kernel(length_scale = 1)
  mk <- batch_kernel(k, batch_size = c(K = 3, O = 2), batch_axes = list(variance = "K"))
  result <- evaluate(mk, c(0, 1))

  for (kk in 1:3) {
    reference <- evaluate(variance_kernel(variance = c(10, 20, 30)[kk]) * se_kernel(length_scale = 1), c(0, 1))
    expect_equal(result[, , kk, 1], reference, tolerance = 1e-12)
    expect_equal(result[, , kk, 2], reference, tolerance = 1e-12)
  }
})

test_that("multi-axis batch_kernel leaves a hyperparameter absent from batch_axes fully shared", {
  k <- variance_kernel(variance = 2) * se_kernel(length_scale = 1)
  mk <- batch_kernel(k, batch_size = c(K = 2, O = 2), batch_axes = list())
  result <- evaluate(mk, c(0, 1))
  reference <- evaluate(k, c(0, 1))
  for (kk in 1:2) for (o in 1:2) {
    expect_equal(result[, , kk, o], reference, tolerance = 1e-12)
  }
})

test_that("multi-axis batch_kernel requires an unnamed batch_size to stay scalar", {
  expect_error(
    batch_kernel(se_kernel(length_scale = 1), batch_size = c(3, 2)),
    "named"
  )
})

test_that("multi-axis batch_kernel requires explicit batch_axes (no auto-detection)", {
  k <- variance_kernel(variance = 1:6) * se_kernel(length_scale = 1)
  expect_error(
    batch_kernel(k, batch_size = c(K = 3, O = 2)),
    "explicit"
  )
})

test_that("multi-axis batch_kernel rejects an unknown axis name in batch_axes", {
  k <- variance_kernel(variance = 1:6) * se_kernel(length_scale = 1)
  expect_error(
    batch_kernel(k, batch_size = c(K = 3, O = 2), batch_axes = list(variance = c("K", "Z"))),
    "unknown axis"
  )
})

test_that("multi-axis batch_kernel rejects a hyperparameter length mismatched with its declared axes", {
  k <- variance_kernel(variance = 1:5) * se_kernel(length_scale = 1) # needs 6, has 5
  expect_error(
    batch_kernel(k, batch_size = c(K = 3, O = 2), batch_axes = list(variance = c("K", "O"))),
    "exactly 6"
  )
})

test_that("multi-axis batch_kernel re-validates its axis contract after kupdate()", {
  k <- variance_kernel(variance = 1:6) * se_kernel(length_scale = 1)
  mk <- batch_kernel(k, batch_size = c(K = 3, O = 2),
                      batch_axes = list(variance = c("K", "O")))
  bad <- kupdate(mk, variance = 1:4) # no longer 6 values
  expect_error(evaluate(bad, c(0, 1)), "exactly 6")
})

test_that("multi-axis batch_kernel with 3 axes matches MIMOSA's (T, K, O) task-kernel shape", {
  # Mirrors build_task_kernel()'s worst case (shared_task_hps = FALSE,
  # cluster_hps_in_tasks = TRUE, shared_output_hps = FALSE): task-kernel
  # hyperparameters vary independently per (task, cluster, output).
  n_task <- 2
  n_cluster <- 2
  n_output <- 3
  n_total <- n_task * n_cluster * n_output
  k <- variance_kernel(variance = seq_len(n_total)) * se_kernel(length_scale = 1)
  mk <- batch_kernel(k, batch_size = c(T = n_task, K = n_cluster, O = n_output),
                      batch_axes = list(variance = c("T", "K", "O")))
  result <- evaluate(mk, c(0, 1))
  expect_equal(dim(result), c(2, 2, n_task, n_cluster, n_output))

  # Spot check flat index for (T=2, K=1, O=3): 2 + (1-1)*2 + (3-1)*2*2 = 10
  expected <- evaluate(variance_kernel(variance = 10) * se_kernel(length_scale = 1), c(0, 1))
  expect_equal(result[, , 2, 1, 3], expected, tolerance = 1e-12)
})

# Regression tests for 5 correctness bugs found in code review of the
# multi-axis contract (dev/ameliorations/design/07): an empty per-hyperparameter axis vector
# used to silently freeze that hyperparameter to its first element for
# every combination; a repeated axis name within one hyperparameter's own
# axis vector used to double-weight that axis and read out-of-range,
# silently NA, indices; a duplicate hyperparameter key across two
# `batch_axes` list entries used to silently keep only the first (`[[`
# lookup semantics) with no error; a single-input-point (1x1)
# self-covariance used to lose its array shape entirely via
# `simplify2array()`'s length-1 collapse; and `block_diag_kernel()` had no
# guard against a multi-axis `nb_blocks`, crashing in `evaluate()` instead
# of failing at construction with a clear message.

test_that("multi-axis batch_kernel rejects an empty per-hyperparameter axis vector", {
  k <- variance_kernel(variance = 1:6) * se_kernel(length_scale = 1)
  expect_error(
    batch_kernel(k, batch_size = c(K = 3, O = 2), batch_axes = list(variance = character(0))),
    "empty"
  )
})

test_that("multi-axis batch_kernel dedupes a repeated axis name within one hyperparameter's vector", {
  k <- variance_kernel(variance = 1:6) * se_kernel(length_scale = 1)
  mk <- batch_kernel(k, batch_size = c(K = 3, O = 2), batch_axes = list(variance = c("K", "O", "K")))
  result <- evaluate(mk, c(0, 1))
  expect_false(anyNA(result))
  # combo (K=1, O=2) -> flat index 1 + (1-1) + (2-1)*3 = 4
  expected <- evaluate(variance_kernel(variance = 4) * se_kernel(length_scale = 1), c(0, 1))
  expect_equal(result[, , 1, 2], expected, tolerance = 1e-12)
})

test_that("multi-axis batch_kernel rejects a duplicate hyperparameter key in batch_axes", {
  k <- variance_kernel(variance = 1:6) * se_kernel(length_scale = 1)
  expect_error(
    batch_kernel(k, batch_size = c(K = 3, O = 2),
                 batch_axes = list(variance = c("K", "O"), variance = c("K"))),
    "more than one entry"
  )
})

test_that("multi-axis batch_kernel keeps the (1, 1, ...) shape for a single input point", {
  k <- variance_kernel(variance = 1:6) * se_kernel(length_scale = 1)
  mk <- batch_kernel(k, batch_size = c(K = 3, O = 2), batch_axes = list(variance = c("K", "O")))
  result <- evaluate(mk, 0)
  expect_equal(dim(result), c(1, 1, 3, 2))
  # combo (K=2, O=1) -> flat index 2; combo (K=1, O=2) -> flat index 4
  expect_equal(result[1, 1, 2, 1], 2, tolerance = 1e-12)
  expect_equal(result[1, 1, 1, 2], 4, tolerance = 1e-12)
})

test_that("block_diag_kernel rejects a non-scalar (multi-axis) nb_blocks", {
  k <- variance_kernel(variance = 1:6) * se_kernel(length_scale = 1)
  expect_error(
    block_diag_kernel(k, nb_blocks = c(K = 3, O = 2), batch_axes = list(variance = c("K", "O"))),
    "single positive integer"
  )
})

test_that("input_specific_param_kernel validates hyperparameter lengths against input_size", {
  expect_error(
    input_specific_param_kernel(white_noise_kernel(noise = c(0.1, 0.5)), input_size = 3),
    "input_size"
  )
})

# The self-covariance contract: whether the requested matrix is a
# self-covariance is decided once, at the level of the user's original
# call, and forwarded explicitly through every recursive evaluate() --
# never re-derived from transformed inputs (a 1-row slice of x1 never
# compares equal to the full x2).

test_that("input_specific_param_kernel keeps na_action = 'mask' invertible", {
  # The masked point's diagonal entry must receive mask_variance -- not 0,
  # which used to leave a whole zero row/column (a singular matrix,
  # breaking the very promise of 'mask').
  k <- input_specific_param_kernel(white_noise_kernel(noise = c(0.1, 0.5, 0.2)),
                                    input_size = 3)
  result <- evaluate(k, c(0, NA, 2), na_action = "mask", mask_variance = 1)
  expect_equal(result[2, 2], 1)
  expect_equal(result[2, c(1, 3)], c(0, 0))
  expect_equal(result[c(1, 3), 2], c(0, 0))
  expect_no_error(chol(result))
})

test_that("input_specific_param_kernel + mask matches the plain kernel for a shared hyperparameter", {
  x <- c(0, NA, 2)
  expect_equal(
    evaluate(input_specific_param_kernel(se_kernel(length_scale = 1), input_size = 3),
             x, na_action = "mask"),
    evaluate(se_kernel(length_scale = 1), x, na_action = "mask"),
    tolerance = 1e-12
  )
})

test_that("batch_kernel forwards the self-covariance status to its slices under na_action = 'mask'", {
  k <- batch_kernel(se_kernel(length_scale = c(1, 2)), batch_size = 2)
  result <- evaluate(k, c(0, NA, 2), na_action = "mask")
  for (b in 1:2) {
    expect_equal(result[2, 2, b], 1) # default mask_variance
    expect_equal(result[2, c(1, 3), b], c(0, 0))
  }
})

test_that("block_kernel's na mask puts mask_variance only in diagonal blocks, 0 in cross blocks", {
  inner <- feature_kernel(length_scale = c(1, 3), length_scale_u = 1, variance = c(2, 1))
  k <- block_kernel(inner, nb_blocks = 2)
  result <- evaluate(k, c(0, NA, 2), na_action = "mask", mask_variance = 7)
  # Point 2 is missing; positions 2 (block 1) and 5 (block 2) are its rows.
  # Diagonal blocks are self-covariances: their own diagonal gets
  # mask_variance...
  expect_equal(result[2, 2], 7)
  expect_equal(result[5, 5], 7)
  # ...but a cross block covaries two *different* processes: the masked
  # point must contribute no cross-information (0), not a self-variance.
  expect_equal(result[2, 5], 0)
  expect_equal(result[5, 2], 0)
})
