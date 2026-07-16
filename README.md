# keRnel

`keRnel` provides composable covariance kernels for Gaussian Process models
in R. Kernels are immutable S3 objects, built around R-native idioms (S3
classes, immutability via copy-on-modify) and combined with arithmetic
operators (`+`, `*`) or wrapped to add batching, automatic relevance
determination, block structure, or input dimension selection.

This package is under active development. The design rationale and the full
inventory of kernels are documented in `dev/ameliorations/design/`:

- `dev/ameliorations/design/01-architecture.Rmd` — full architecture (class system,
  hyperparameter representation, the 4-layer computation model, freeze/
  share, `kupdate()`, validation, package conventions).
- `dev/ameliorations/design/03-analyse-risques.Rmd` — anticipated methodological risks and
  the recommended implementation order.
- `dev/ameliorations/design/06-choix-methodologiques-kernax-vs-keRnel.Rmd` — the
  methodological choices carried over from `kernax` (the Python/JAX sibling
  library this package was ported from), and where the two implementations
  deliberately diverge.

## Status

Following the implementation order from `03-analyse-risques.Rmd`:

- [x] `wrapped_param`, parametrisations, distance functions
- [x] `se_kernel` end-to-end (template for all base kernels)
- [x] Vectorised `distance_matrix()`/`pairwise_matrix()` for stationary
      kernels (brought forward after measuring the naive `map_pairs()` loop
      at scale -- see `dev/ameliorations/design/03-analyse-risques.Rmd`, risk #4)
- [x] Remaining base and dot-product kernels (`matern12_kernel`,
      `matern32_kernel`, `matern52_kernel`, `periodic_kernel`,
      `rational_quadratic_kernel`, `white_noise_kernel`, `feature_kernel`,
      `constant_kernel`, `variance_kernel`, `linear_kernel`,
      `affine_kernel`, `polynomial_kernel`, `sigmoid_kernel`)
- [x] `diagonal_engine` (a fast path for self-covariance diagonal
      computation; see `dev/ameliorations/design/01-architecture.Rmd`), with a vectorised
      `pairwise_diag()`/`diag_distance()` path mirroring
      `pairwise_matrix()`/`distance_matrix()`. Measured: 10,000 points in
      0.02s (`O(n)`) vs. `dense_engine()`'s 0.54s for just 2,000 points
      (`O(n^2)`)
- [x] Operators (`sum_kernel`/`+`, `product_kernel`/`*`) and simple wrappers
      (`exp_kernel`, `log_kernel`), plus the `kernel_interactions()`
      ergonomic helper. `neg_kernel`/`-`/unary negation deliberately *not*
      implemented (explicit user decision -- subtraction of covariance
      kernels is rarely meaningful and was judged unnecessary).
      `dim_additive_kernel()` deferred to the `active_dims_kernel()` step.
- [x] `evaluate()` NA handling (`na_action = c("propagate", "mask")`,
      `mask_variance`) and the memory-size warning guard
      (`memory_warning_threshold`), forwarded through composite/wrapper
      kernels via `...`
- [x] `active_dims_kernel`, `ard_kernel` (the latter takes a kernel
      *constructor*, not a pre-built instance -- builds and `freeze()`s the
      inner kernel itself, so the "the inner length_scale must be 1 and
      non-trainable" usage contract cannot be violated by the caller)
- [x] `batch_kernel`, `input_specific_param_kernel`, `block_kernel`,
      `block_diag_kernel` -- which hyperparameters vary along a wrapper's
      batch/block axis is an explicit, validated contract (auto-detected by
      default, or declared by name via `batch_axes`/`block_axes`/
      `param_axes`), re-checked at every `evaluate()` so a later
      `kupdate()` cannot silently desynchronise it. Nested
      `batch_kernel(batch_kernel(...))`, each level owning its own axis
      (including different sizes per level), is covered by tests. The
      self-covariance status of a call (needed by `na_action = "mask"`) is
      decided once, at the level of the user's call, and threaded
      explicitly through every internal `evaluate()` -- never re-derived
      from an already-transformed slice.
- [x] `sample_hps()` (a fourth consumer of the `flatten_params()`-style
      recursive walk, alongside `get_params()`/`kupdate()`/`freeze()`),
      and four comprehensive vignettes covering every documented behaviour:
  - `vignettes/a-beginner.Rmd`: every base kernel, printing,
    multi-dimensional inputs, `+`/`*` composition, default NA handling.
  - `vignettes/b-intermediate.Rmd`: `get_params()`/`kupdate()`/
    `freeze()`/`unfreeze()`, all 5 parametrisations, the full distance
    catalogue, `exp_kernel()`/`log_kernel()`, `kernel_interactions()`,
    `active_dims_kernel()`, `ard_kernel()`, `sample_hps()`,
    `na_action = "mask"`, writing a custom kernel.
  - `vignettes/c-expert.Rmd`: the 4-layer architecture, `dense_engine()`/
    `diagonal_engine()` (with a real performance measurement), the
    vectorisation internals (`shape()`), `batch_kernel()` (incl. nested
    batching), `input_specific_param_kernel()`,
    `block_kernel()`/`block_diag_kernel()`, the memory warning guard.
  - `vignettes/d-batching.Rmd`: batching as a Gaussian Process modelling
    tool -- GP prior draws under batched hyperparameters, batched inputs,
    the explicit axis contract, nested batching grids, heteroscedastic
    regression via `input_specific_param_kernel()`, missing data across
    structured kernels, and multi-output GPs via `block_kernel()`.

  Writing the vignettes surfaced a real bug (not caught by any existing
  `testthat` test): `pairwise_matrix.stationary_kernel()`/
  `pairwise_diag.stationary_kernel()` called `shape()` unconditionally,
  crashing for a custom kernel that only defines `pairwise()` (a
  perfectly valid kernel, just without the vectorised fast path). Fixed
  with `has_s3_method()`, falling back to the always-correct
  `map_pairs()`/`map_diag()` path when `shape()` is absent. Also led to
  exporting `resolve_distance_function()`, `validate_kernel()`, and
  `shape()`, since the vignettes' "Writing your own kernel" section
  established these as legitimate public extension API, not just
  implementation details.
- [x] `batch_kernel()` multi-axis mode: a single wrapper can now make one
      hyperparameter vary along the full Cartesian product of 2+ named
      axes at once (`batch_size = c(K = 3, O = 2)`, `batch_axes` a named
      list), something nested single-axis `batch_kernel()` calls cannot
      express. Motivated by analysing MIMOSA (the Python/JAX multi-task
      GP framework `keRnel`'s kernels are meant to eventually support an
      R port of) -- see `objectifs/01-hp-sharing-matrix.md` for the
      scenario matrix this unblocks and `dev/ameliorations/design/07-contrat-multi-axes.Rmd`
      for the design rationale.

## Development

```r
devtools::load_all()
devtools::test()
devtools::document()  # regenerates NAMESPACE/man from roxygen comments
```

Benchmarks (not part of the installed package, see
`dev/benchmark/README.md`):

```r
bench::mark(...)  # see dev/benchmark/ for ready-made scenarios
```
