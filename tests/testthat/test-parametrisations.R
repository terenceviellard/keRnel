test_that("identity_parametrisation leaves values unchanged", {
  p <- identity_parametrisation()
  expect_equal(wrap(p, 5), 5)
  expect_equal(unwrap(p, -3), -3)
})

test_that("log_exp_parametrisation wraps/unwraps positive values correctly", {
  p <- log_exp_parametrisation()
  expect_equal(wrap(p, 1), 0)
  expect_equal(unwrap(p, 0), 1)
  expect_equal(unwrap(p, wrap(p, 4.2)), 4.2, tolerance = 1e-8)
})

test_that("log_exp_parametrisation rejects non-positive values", {
  p <- log_exp_parametrisation()
  expect_error(wrap(p, 0), "strictly positive")
  expect_error(wrap(p, -1), "strictly positive")
})

test_that("softplus_parametrisation round-trips positive values", {
  p <- softplus_parametrisation()
  expect_equal(unwrap(p, wrap(p, 3.3)), 3.3, tolerance = 1e-8)
  expect_error(wrap(p, -1), "strictly positive")
})

test_that("softplus and log_exp parametrisations agree at the round trip but not on the wrapped value", {
  x <- 2.7
  expect_equal(
    unwrap(softplus_parametrisation(), wrap(softplus_parametrisation(), x)),
    unwrap(log_exp_parametrisation(), wrap(log_exp_parametrisation(), x)),
    tolerance = 1e-8
  )
  expect_false(isTRUE(all.equal(
    wrap(softplus_parametrisation(), x),
    wrap(log_exp_parametrisation(), x)
  )))
})

test_that("bounded_parametrisation round-trips values strictly within bounds", {
  p <- bounded_parametrisation(lower = 0, upper = 10)
  expect_equal(unwrap(p, wrap(p, 4)), 4, tolerance = 1e-8)
})

test_that("bounded_parametrisation rejects values outside (lower, upper)", {
  p <- bounded_parametrisation(lower = 0, upper = 10)
  expect_error(wrap(p, -1), "strictly within")
  expect_error(wrap(p, 0), "strictly within")
  expect_error(wrap(p, 10), "strictly within")
})

test_that("bounded_parametrisation requires lower < upper", {
  expect_error(bounded_parametrisation(5, 5), "strictly less")
  expect_error(bounded_parametrisation(5, 2), "strictly less")
})

test_that("parametrisation_chain applies wrap() in order and unwrap() in reverse order", {
  chain <- parametrisation_chain(log_exp_parametrisation(), identity_parametrisation())
  direct <- log_exp_parametrisation()
  expect_equal(wrap(chain, 2), wrap(direct, 2))
  expect_equal(unwrap(chain, wrap(chain, 2)), 2, tolerance = 1e-8)
})

test_that("parametrisation_chain with no parametrisations behaves like identity", {
  chain <- parametrisation_chain()
  expect_equal(wrap(chain, 5), 5)
  expect_equal(unwrap(chain, 5), 5)
})
