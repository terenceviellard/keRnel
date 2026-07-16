# Meta-test: the actual safety net behind validate_kernel() (see
# dev/ameliorations/design/01-architecture.Rmd, "Garde-fous structurels sans S4"). Every
# constructor of a `stationary_kernel`/`dot_product_kernel` must be listed
# here so the contract is checked systematically rather than trusted by
# convention alone.
#
# TODO: extend these lists as new stationary/dot-product kernels are added.

stationary_constructors <- list(
  se_kernel = function() se_kernel(length_scale = 1),
  rbf_kernel = function() rbf_kernel(length_scale = 1),
  matern12_kernel = function() matern12_kernel(length_scale = 1),
  matern32_kernel = function() matern32_kernel(length_scale = 1),
  matern52_kernel = function() matern52_kernel(length_scale = 1),
  periodic_kernel = function() periodic_kernel(length_scale = 1, period = 1),
  rational_quadratic_kernel = function() rational_quadratic_kernel(length_scale = 1, alpha = 1),
  white_noise_kernel = function() white_noise_kernel(noise = 1),
  feature_kernel = function() feature_kernel(length_scale = 1, length_scale_u = 1, variance = 1)
)

dot_product_constructors <- list(
  linear_kernel = function() linear_kernel(slope_var = 1),
  affine_kernel = function() affine_kernel(slope_var = 1, offset = 0),
  polynomial_kernel = function() polynomial_kernel(degree = 2, gamma = 1, constant = 0),
  sigmoid_kernel = function() sigmoid_kernel(alpha = 1, constant = 0)
)

test_that("every stationary_kernel constructor produces a function-valued distance_function", {
  for (name in names(stationary_constructors)) {
    kernel <- stationary_constructors[[name]]()
    expect_true(
      is.function(kernel$distance_function),
      info = paste(name, "must have a function-valued distance_function")
    )
  }
})

test_that("every dot_product_kernel constructor produces a function-valued distance_function", {
  for (name in names(dot_product_constructors)) {
    kernel <- dot_product_constructors[[name]]()
    expect_true(
      is.function(kernel$distance_function),
      info = paste(name, "must have a function-valued distance_function")
    )
  }
})

test_that("validate_kernel() rejects a hand-built stationary_kernel missing distance_function", {
  broken <- structure(
    list(length_scale = wrapped_param(1)),
    class = c("broken_kernel", "stationary_kernel", "kernel")
  )
  expect_error(validate_kernel(broken), "distance_function")
})

test_that("validate_kernel() rejects a hand-built dot_product_kernel missing distance_function", {
  broken <- structure(
    list(slope_var = wrapped_param(1)),
    class = c("broken_kernel", "dot_product_kernel", "kernel")
  )
  expect_error(validate_kernel(broken), "distance_function")
})

operator_constructors <- list(
  sum_kernel = function() sum_kernel(se_kernel(length_scale = 1), se_kernel(length_scale = 1)),
  product_kernel = function() product_kernel(se_kernel(length_scale = 1), se_kernel(length_scale = 1))
)

wrapper_constructors <- list(
  exp_kernel = function() exp_kernel(se_kernel(length_scale = 1)),
  log_kernel = function() log_kernel(se_kernel(length_scale = 1)),
  active_dims_kernel = function() active_dims_kernel(se_kernel(length_scale = 1), active_dims = 1),
  ard_kernel = function() ard_kernel(se_kernel, length_scales = c(1, 2)),
  batch_kernel = function() batch_kernel(se_kernel(length_scale = c(1, 2)), batch_size = 2),
  input_specific_param_kernel = function() {
    input_specific_param_kernel(se_kernel(length_scale = 1), input_size = 2)
  },
  block_diag_kernel = function() block_diag_kernel(se_kernel(length_scale = c(1, 2)), nb_blocks = 2),
  block_kernel = function() block_kernel(se_kernel(length_scale = 1), nb_blocks = 2)
)

test_that("every operator_kernel constructor produces exactly 2 kernel-valued subkernels", {
  for (name in names(operator_constructors)) {
    kernel <- operator_constructors[[name]]()
    expect_length(kernel$subkernels, 2)
    expect_true(all(vapply(kernel$subkernels, inherits, logical(1), what = "kernel")))
  }
})

test_that("every wrapper_kernel constructor produces exactly 1 kernel-valued subkernel", {
  for (name in names(wrapper_constructors)) {
    kernel <- wrapper_constructors[[name]]()
    expect_length(kernel$subkernels, 1)
    expect_true(inherits(kernel$subkernels[[1]], "kernel"))
  }
})

# Generic behavioural contract (finding #7): every base kernel must
# implement pairwise() (or override evaluate() itself), and have an
# `engine`; every operator_kernel/wrapper_kernel must override evaluate().
# Checked once, generically, via S3 method presence -- not per concrete
# class -- so a new kernel is covered automatically.

test_that("validate_kernel() rejects a base kernel with neither pairwise() nor evaluate()", {
  broken <- structure(list(engine = dense_engine()), class = c("broken_kernel", "kernel"))
  expect_error(validate_kernel(broken), "pairwise")
})

test_that("validate_kernel() rejects a base kernel with no `engine`", {
  broken <- structure(list(), class = c("constant_kernel", "kernel"))
  broken$value <- wrapped_param(1, identity_parametrisation())
  expect_error(validate_kernel(broken), "engine")
})

test_that("validate_kernel() rejects an operator_kernel that does not override evaluate()", {
  broken <- structure(
    list(subkernels = list(se_kernel(length_scale = 1), se_kernel(length_scale = 1))),
    class = c("broken_op_kernel", "operator_kernel", "kernel")
  )
  expect_error(validate_kernel(broken), "override `evaluate")
})

test_that("validate_kernel() rejects a wrapper_kernel that does not override evaluate()", {
  broken <- structure(
    list(subkernels = list(se_kernel(length_scale = 1))),
    class = c("broken_wrap_kernel", "wrapper_kernel", "kernel")
  )
  expect_error(validate_kernel(broken), "override `evaluate")
})

test_that("every base kernel constructor passes the generic pairwise()/engine contract", {
  base_constructors <- c(stationary_constructors, dot_product_constructors, list(
    constant_kernel = function() constant_kernel(value = 1),
    variance_kernel = function() variance_kernel(variance = 1)
  ))
  for (name in names(base_constructors)) {
    expect_no_error(validate_kernel(base_constructors[[name]]()))
  }
})

test_that("every operator_kernel/wrapper_kernel constructor passes the generic evaluate() contract", {
  for (name in names(operator_constructors)) {
    expect_no_error(validate_kernel(operator_constructors[[name]]()))
  }
  for (name in names(wrapper_constructors)) {
    expect_no_error(validate_kernel(wrapper_constructors[[name]]()))
  }
})
