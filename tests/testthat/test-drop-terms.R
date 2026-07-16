test_that("drop_terms() removes a matching term from a sum", {
  k <- se_kernel(length_scale = 1) + white_noise_kernel(noise = 0.1)
  result <- drop_terms(k, "white_noise_kernel")
  expect_equal(result, se_kernel(length_scale = 1))
})

test_that("drop_terms() removes only the matching term(s), keeping the rest of the sum", {
  k <- se_kernel(length_scale = 1) + variance_kernel(2) + white_noise_kernel(noise = 0.1)
  result <- drop_terms(k, "white_noise_kernel")
  expect_equal(evaluate(result, c(0, 1, 2)), evaluate(se_kernel(length_scale = 1) + variance_kernel(2), c(0, 1, 2)))
})

test_that("drop_terms() accepts a predicate function instead of a class name", {
  k <- se_kernel(length_scale = 1) + white_noise_kernel(noise = 0.1)
  result <- drop_terms(k, function(x) inherits(x, "white_noise_kernel"))
  expect_equal(result, se_kernel(length_scale = 1))
})

test_that("drop_terms() errors when the match is a factor of a product", {
  k <- variance_kernel(2) * se_kernel(length_scale = 1)
  expect_error(drop_terms(k, "variance_kernel"), "factor")
})

test_that("drop_terms() errors when every term matches (empty kernel)", {
  k <- se_kernel(length_scale = 1)
  expect_error(drop_terms(k, "se_kernel"), "empty kernel")
})

test_that("drop_terms() errors when the match is inside an unsupported wrapper", {
  k <- exp_kernel(white_noise_kernel(noise = 0.1))
  expect_error(drop_terms(k, "white_noise_kernel"), "ambiguous")
})
