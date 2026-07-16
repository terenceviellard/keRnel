test_that("precompute()/eval_plan() match evaluate() for a simple stationary kernel", {
  k <- se_kernel(length_scale = 1)
  plan <- precompute(k, c(0, 1, 2))
  expect_equal(eval_plan(plan, se_kernel(length_scale = 2)), evaluate(se_kernel(length_scale = 2), c(0, 1, 2)))
})

test_that("eval_plan() errors clearly when kernel's structure does not match the plan", {
  k3 <- se_kernel(length_scale = 1) + white_noise_kernel(noise = 0.1)
  plan3 <- precompute(k3, c(0, 1, 2))
  k4 <- se_kernel(length_scale = 1) * white_noise_kernel(noise = 0.1) # same leaves, different operator
  expect_error(eval_plan(plan3, k4), "does not match the tree structure")
})

test_that("eval_plan() errors clearly when a leaf's class differs from the plan's", {
  k1 <- se_kernel(length_scale = 1)
  plan1 <- precompute(k1, c(0, 1, 2))
  expect_error(eval_plan(plan1, periodic_kernel(length_scale = 1, period = 2)), "does not match the tree structure")
})

test_that("a fully-fallback tree (nothing cacheable) still matches evaluate()", {
  k <- affine_kernel(slope_var = 2, offset = 0.5)
  plan <- precompute(k, c(0, 1, 2))
  k2 <- affine_kernel(slope_var = 3, offset = 0.5)
  expect_equal(eval_plan(plan, k2), evaluate(k2, c(0, 1, 2)))
})

test_that("2+ levels of nesting mixing cacheable and fallback subtrees still match evaluate()", {
  inner <- (se_kernel(length_scale = 1) + affine_kernel(slope_var = 1, offset = 0)) +
    (white_noise_kernel(noise = 0.1) + affine_kernel(slope_var = 2, offset = 1))
  plan <- precompute(inner, c(0, 1, 2, 3))
  inner2 <- (se_kernel(length_scale = 2) + affine_kernel(slope_var = 1, offset = 0)) +
    (white_noise_kernel(noise = 0.2) + affine_kernel(slope_var = 2, offset = 1))
  expect_equal(eval_plan(plan, inner2), evaluate(inner2, c(0, 1, 2, 3)))
})

test_that("eval_plan() reproduces NA handling exactly, including through a nested fallback subtree", {
  inner <- (se_kernel(length_scale = 1) + affine_kernel(slope_var = 1, offset = 0))
  x <- c(0, NA, 2, 3)
  plan <- precompute(inner, x)
  for (na_mode in c("propagate", "mask")) {
    expect_equal(eval_plan(plan, inner, na_action = na_mode), evaluate(inner, x, na_action = na_mode))
  }
})

test_that("eval_plan()'s operator combination also rejects a dense/diagonal engine mismatch", {
  inner <- se_kernel(length_scale = c(1, 2)) + se_kernel(length_scale = 1, engine = diagonal_engine())
  k <- batch_kernel(inner, batch_size = 2)
  expect_error(evaluate(k, c(0, 1, 2), c(0, 5, 100)), "dense_engine\\(\\).*diagonal_engine\\(\\)")
})

test_that("get_free_params()/set_free_params() work correctly on a batch_kernel using the plan-based path", {
  k <- batch_kernel(se_kernel(length_scale = c(1, 2, 3)), batch_size = 3)
  free <- get_free_params(k)
  expect_length(free, 3)
  k2 <- set_free_params(k, free + 1)
  for (b in 1:3) {
    expect_equal(get_params(slice_kernel(k2$subkernels[[1]], b, "length_scale"))[["length_scale"]],
                 exp(log(b) + 1), tolerance = 1e-8)
  }
})
