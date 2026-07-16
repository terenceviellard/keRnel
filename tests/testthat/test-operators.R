# sum_kernel: k(x1, x2) = k1(x1, x2) + k2(x1, x2)

test_that("sum_kernel adds the pointwise results of its two subkernels", {
  k1 <- se_kernel(length_scale = 1)
  k2 <- white_noise_kernel(noise = 0.5)
  k <- sum_kernel(k1, k2)
  expect_equal(pairwise(k1, c(0), c(0)) + pairwise(k2, c(0), c(0)),
               evaluate(k, c(0))[1, 1])
})

test_that("`+` builds a sum_kernel", {
  k1 <- se_kernel(length_scale = 1)
  k2 <- linear_kernel(slope_var = 1)
  k <- k1 + k2
  expect_true(inherits(k, "sum_kernel"))
  expect_equal(evaluate(k, c(0, 1, 2)), evaluate(k1, c(0, 1, 2)) + evaluate(k2, c(0, 1, 2)))
})

test_that("`+` coerces a bare number to a constant_kernel", {
  k1 <- se_kernel(length_scale = 1)
  k <- k1 + 2
  expect_true(inherits(k$subkernels[[2]], "constant_kernel"))
  expect_equal(evaluate(k, c(0, 1)), evaluate(k1, c(0, 1)) + 2)
})

test_that("`+` is commutative in operand coercion (number + kernel)", {
  k1 <- se_kernel(length_scale = 1)
  k <- 2 + k1
  expect_true(inherits(k, "sum_kernel"))
  expect_equal(evaluate(k, c(0, 1)), 2 + evaluate(k1, c(0, 1)))
})

test_that("sum_kernel's class vector reflects its category hierarchy", {
  k <- se_kernel(length_scale = 1) + linear_kernel(slope_var = 1)
  expect_equal(class(k), c("sum_kernel", "operator_kernel", "kernel"))
})

test_that("sum_kernel errors clearly when combining a dense_engine sub-kernel with a diagonal_engine one", {
  k <- se_kernel(length_scale = 1) + se_kernel(length_scale = 1, engine = diagonal_engine())
  expect_error(evaluate(k, c(0, 1, 2), c(0, 5, 100)), "dense_engine\\(\\).*diagonal_engine\\(\\)")
})

test_that("sum_kernel of two diagonal_engine sub-kernels combines correctly (both vectors)", {
  k <- se_kernel(length_scale = 1, engine = diagonal_engine()) +
    white_noise_kernel(noise = 0.1, engine = diagonal_engine())
  result <- evaluate(k, c(0, 1, 2), c(0, 5, 100))
  expect_false(is.matrix(result))
  expect_length(result, 3)
})

# product_kernel: k(x1, x2) = k1(x1, x2) * k2(x1, x2)

test_that("product_kernel multiplies the pointwise results of its two subkernels", {
  k1 <- variance_kernel(variance = 2)
  k2 <- se_kernel(length_scale = 1)
  k <- product_kernel(k1, k2)
  expect_equal(evaluate(k, c(0, 1)), 2 * evaluate(k2, c(0, 1)))
})

test_that("product_kernel errors clearly when combining a dense_engine sub-kernel with a diagonal_engine one", {
  k <- se_kernel(length_scale = 1) * se_kernel(length_scale = 1, engine = diagonal_engine())
  expect_error(evaluate(k, c(0, 1, 2), c(0, 5, 100)), "dense_engine\\(\\).*diagonal_engine\\(\\)")
})

test_that("`*` builds a product_kernel", {
  k1 <- variance_kernel(variance = 3)
  k2 <- se_kernel(length_scale = 1)
  k <- k1 * k2
  expect_true(inherits(k, "product_kernel"))
  expect_equal(evaluate(k, c(0, 1, 2)), evaluate(k1, c(0, 1, 2)) * evaluate(k2, c(0, 1, 2)))
})

test_that("`*` coerces a bare number to a constant_kernel", {
  k1 <- se_kernel(length_scale = 1)
  k <- k1 * 2
  expect_true(inherits(k$subkernels[[2]], "constant_kernel"))
  expect_equal(evaluate(k, c(0, 1)), evaluate(k1, c(0, 1)) * 2)
})

# Operators recurse through get_params()/kupdate()/freeze() for free, via
# the uniform $subkernels storage -- no new code was needed for these.

test_that("get_params() collects hyperparameters from both subkernels", {
  k <- se_kernel(length_scale = 2) + white_noise_kernel(noise = 0.3)
  params <- get_params(k)
  expect_equal(unname(params), c(2, 0.3), tolerance = 1e-10)
})

test_that("kupdate() updates a hyperparameter wherever it appears in the tree", {
  k <- se_kernel(length_scale = 1) * se_kernel(length_scale = 1)
  k2 <- kupdate(k, length_scale = 5)
  expect_equal(unname(get_params(k2)), c(5, 5))
})

test_that("freeze() excludes a subkernel's hyperparameter from get_trainable_params()", {
  k <- se_kernel(length_scale = 1) + white_noise_kernel(noise = 0.5)
  k <- freeze(k, "noise")
  expect_length(get_trainable_params(k), 1)
  expect_length(get_params(k), 2)
})

# Unsupported operators

test_that("`-` is explicitly not supported for kernels", {
  k1 <- se_kernel(length_scale = 1)
  k2 <- se_kernel(length_scale = 2)
  expect_error(k1 - k2, "not defined")
})

test_that("unary `-` is explicitly not supported for kernels", {
  k1 <- se_kernel(length_scale = 1)
  expect_error(-k1, "not defined")
})

# Validation

test_that("validate_kernel() rejects a hand-built operator_kernel with a non-kernel subkernel", {
  broken <- structure(
    list(subkernels = list(se_kernel(length_scale = 1), 4)),
    class = c("broken_kernel", "operator_kernel", "kernel")
  )
  expect_error(validate_kernel(broken), "must itself be a `kernel`")
})

test_that("validate_kernel() rejects a hand-built operator_kernel with the wrong number of subkernels", {
  broken <- structure(
    list(subkernels = list(se_kernel(length_scale = 1))),
    class = c("broken_kernel", "operator_kernel", "kernel")
  )
  expect_error(validate_kernel(broken), "exactly 2")
})

# Printing

test_that("a sum_kernel's format() composes the format() of its subkernels", {
  k <- se_kernel(length_scale = 1) + linear_kernel(slope_var = 2)
  expect_match(format(k), "\\+")
})

test_that("a product_kernel parenthesises a sum_kernel operand", {
  k1 <- se_kernel(length_scale = 1)
  k2 <- linear_kernel(slope_var = 1)
  k3 <- white_noise_kernel(noise = 1)
  k <- (k1 + k2) * k3
  expect_match(format(k), "^\\(.*\\+.*\\) \\*")
})

# as_kernel() hardening: only a kernel or a single finite numeric value

test_that("`+`/`*` reject NULL, NA, and non-numeric operands instead of silently constructing a broken kernel", {
  k1 <- se_kernel(length_scale = 1)
  expect_error(k1 + NULL, "kernel.*finite numeric")
  expect_error(k1 + NA, "kernel.*finite numeric")
  expect_error(k1 + "x", "kernel.*finite numeric")
  expect_error(k1 * NULL, "kernel.*finite numeric")
  expect_error(k1 + c(1, 2), "kernel.*finite numeric")
})
