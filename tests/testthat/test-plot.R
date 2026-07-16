test_that("plot_kernel() returns a ggplot object built from the kernel's covariance matrix", {
  testthat::skip_if_not_installed("ggplot2")
  k <- se_kernel(length_scale = 1)
  p <- plot_kernel(k, c(0, 1, 2))
  expect_s3_class(p, "ggplot")
  expect_equal(p$data$k, as.vector(evaluate(k, c(0, 1, 2))))
})
