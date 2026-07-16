test_that("format.kernel shows the kernel class and its hyperparameters", {
  k <- se_kernel(length_scale = 2)
  expect_match(format(k), "^se_kernel\\(")
  expect_match(format(k), "length_scale = 2")
})

test_that("format.kernel reflects a frozen hyperparameter", {
  k <- freeze(se_kernel(length_scale = 2), "length_scale")
  expect_match(format(k), "frozen")
})

test_that("print.kernel prints the formatted kernel and returns it invisibly", {
  k <- se_kernel(length_scale = 2)
  expect_output(print(k), "se_kernel\\(length_scale = 2\\)")
  expect_identical(withVisible(print(k))$visible, FALSE)
})

test_that("kernel_tree() prints one indented line per node of the $subkernels tree", {
  k <- variance_kernel(0.4) * se_kernel(length_scale = 2)
  out <- capture.output(kernel_tree(k))
  expect_length(out, 3) # product_kernel, variance_kernel, se_kernel
  expect_match(out[1], "^product_kernel")
  expect_match(out[2], "^  variance_kernel \\(variance = 0.4\\)")
  expect_match(out[3], "^  se_kernel \\(length_scale = 2\\)")
})

test_that("kernel_tree() shows a frozen hyperparameter and returns the kernel invisibly", {
  k <- freeze(white_noise_kernel(noise = 0.05), "noise")
  expect_output(kernel_tree(k), "\\[frozen\\]")
  expect_identical(withVisible(kernel_tree(k))$visible, FALSE)
})
