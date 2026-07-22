# keRnel 0.1.0

Initial release. `keRnel` provides a collection of composable covariance
kernel functions for Gaussian Process models: immutable S3 objects that can
be combined with arithmetic operators (`+`, `*`) and wrapped to add
batching, automatic relevance determination, block structure, or input
dimension selection.

## Performance

* `batch_kernel()` now evaluates through an internal cached distance plan
  whenever `batch_over_inputs = FALSE`: the distance matrix of any
  stationary/dot-product leaf is computed once and reused across batch
  elements instead of being recomputed per slice, with an automatic,
  per-subtree fallback for any kernel the plan cannot cache. Measured
  ~1.9x faster than the equivalent hand-written loop.
* The underlying mechanism is exposed directly via `precompute()`/
  `eval_plan()`, for callers that re-evaluate the same kernel structure
  many times over a fixed set of points (e.g. an optimiser, or an EM/VEM
  training loop).

## New functions

* `evaluate_named()`: like `evaluate()`, but propagates `x1`/`x2`'s point
  identity (`rownames()`/`names()`) onto the result's dimnames.
* `robust_chol()` / `robust_chol_inv()`: Cholesky factor / inverse with
  adaptive-jitter recovery from a near-singular matrix, dimnames preserved.
* `drop_terms()`: removes a matching term from a kernel sum, rebuilding the
  tree rather than approximating removal numerically.
* `plot_kernel()`: a `ggplot2` heatmap of a kernel's covariance matrix in
  one call.
* `get_free_params()` / `set_free_params()`: read/write a kernel's
  non-frozen hyperparameters directly in the unconstrained (wrapped) space,
  positionally rather than by name -- avoids the non-bit-exact
  `wrap()`/`unwrap()` round trip an optimiser working in unconstrained
  space would otherwise pay every iteration, and correctly handles two
  hyperparameters sharing a bare name at different points in the tree
  (which `kupdate()` cannot set independently).
* `kernel_grad()`: the analytic gradient of `evaluate(kernel, x1, x2)` with
  respect to each non-frozen hyperparameter, for every base kernel, any
  `sum_kernel()`/`product_kernel()` composition of them (including
  `kernel_interactions()`), `exp_kernel()`, `log_kernel()`,
  `active_dims_kernel()`, and `ard_kernel()` -- avoids the `2p + 1`
  finite-difference evaluations per step a gradient-free optimiser needs.
  `shape_grad()`/`shape_grad_d()`/`engine_grad()` are the underlying
  per-kernel/per-engine building blocks, mirroring `shape()`/
  `engine_eval()`; `ard_kernel()`'s gradient additionally requires its
  wrapped kernel's `distance_function` to be one of
  `squared_euclidean_distance()`/`euclidean_distance()`/
  `manhattan_distance()`/`dot_product()`. Not yet supported: `batch_kernel()`,
  `block_kernel()`, `block_diag_kernel()`, `input_specific_param_kernel()`.
* `kernel_grad_free()`: like `kernel_grad()`, chain-ruled (via the new
  `unwrap_grad()`) into the free/unconstrained space `get_free_params()`/
  `set_free_params()` work in -- positional and unnamed, so it stays
  correct even when two sub-kernels share a bare hyperparameter name (a
  case `kupdate()`/`kernel_grad()`'s name-based APIs cannot address
  independently). The intended pairing for an unconstrained optimiser
  (e.g. `optim(method = "L-BFGS-B")` without box constraints, since
  positivity is already guaranteed by each hyperparameter's
  parametrisation).

## Bug fixes

* `evaluate()` now treats `Inf`/`-Inf` in the input like `NA` (subject to
  `na_action`), rather than letting them produce a silent `NaN` in the
  covariance matrix.
* `as_matrix_input()` rejects a matrix with 0 columns with a clear error,
  instead of a misleading "all points identical" result.
* `sum_kernel()`/`product_kernel()` (`+`/`*`) now error clearly when
  combining a sub-kernel using `dense_engine()` (a matrix result) with one
  using `diagonal_engine()` (a vector result), instead of silently
  combining them via R's vector-recycling rules -- a result that was not a
  meaningful covariance and could pass `check_kernel_contract()` by
  accident under self-covariance.
* `eval_plan()` now errors clearly if `kernel`'s tree structure does not
  match the plan `precompute()` built (e.g. an operator changed from `+`
  to `*`) -- previously it silently reused the plan's original structure,
  producing a wrong result with no warning.
* `eval_plan()` called directly (not through `batch_kernel()`) now
  correctly defaults `.self_covariance` from the `x1`/`x2` originally
  passed to `precompute()`, instead of always assuming `FALSE` -- this
  affected `na_action = "mask"` masking when calling `eval_plan()` outside
  of `batch_kernel()`.
* `kupdate()` now refuses to silently collapse a multi-valued (batch/block
  axis) hyperparameter down to a single shared value -- previously
  `kupdate(k, length_scale = 1)` on a 3-element batch axis would silently
  make all 3 batch elements identical, with no error.
* `evaluate.batch_kernel()`/`block_diag_kernel()` no longer lose the array
  shape when evaluating a single input point: `block_diag_kernel(inner,
  nb_blocks = 1)` on one point used to crash outright; `batch_kernel()` on
  one point used to silently return a bare vector instead of a
  `(1, 1, batch_size)` array.

## Documentation

* `?periodic_kernel` documents the numerical aliasing that occurs for an
  extremely small `period`.

## Infrastructure

* Continuous integration (`R CMD check` across OS/R versions, lint, test
  coverage) and a `.lintr` configuration.
