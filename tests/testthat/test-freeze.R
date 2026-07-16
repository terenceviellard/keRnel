test_that("freeze() marks the named parameter as frozen", {
  k <- freeze(se_kernel(length_scale = 1), "length_scale")
  expect_true(k$length_scale$frozen)
})

test_that("unfreeze() reverses freeze()", {
  k <- se_kernel(length_scale = 1)
  k_frozen <- freeze(k, "length_scale")
  k_unfrozen <- unfreeze(k_frozen, "length_scale")
  expect_false(k_unfrozen$length_scale$frozen)
})

test_that("freeze()/unfreeze() error on an unknown hyperparameter name", {
  k <- se_kernel(length_scale = 1)
  expect_error(freeze(k, "not_a_param"), "No such hyperparameter")
  expect_error(unfreeze(k, "not_a_param"), "No such hyperparameter")
})

test_that("freeze() with no names returns the kernel unchanged", {
  k <- se_kernel(length_scale = 1)
  expect_equal(freeze(k), k)
})

test_that("freezing does not affect the constrained value, only optimisability", {
  k <- freeze(se_kernel(length_scale = 3), "length_scale")
  expect_equal(get_params(k)[["length_scale"]], 3, tolerance = 1e-8)
})
