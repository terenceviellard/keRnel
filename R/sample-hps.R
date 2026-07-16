# A fourth consumer of the same recursive tree-walk as flatten_params()/
# kupdate()/freeze() (see dev/ameliorations/design/01-architecture.Rmd) -- randomly
# re-initialising named hyperparameters from a prior distribution, useful
# for multi-start optimisation of a GP fit.

#' Randomly sample a kernel's hyperparameters from a prior distribution
#'
#' @description
#' Returns a copy of `kernel` with the named hyperparameter(s) resampled
#' from `priors`, in the kernel's natural (constrained) space, applied
#' wherever each hyperparameter exists in the (possibly nested) kernel tree
#' -- the same recursive walk as [kupdate()] (via [flatten_params()]).
#'
#' Frozen hyperparameters (see [freeze()]) are skipped, for consistency with
#' [get_trainable_params()]: resampling a value that will then be excluded
#' from optimisation is rarely useful, and risks silently overwriting a
#' value the caller deliberately fixed (e.g. [ard_kernel()]'s inner
#' `length_scale`, always frozen at `1`).
#'
#' `sample_hps()` relies on R's ambient RNG state -- call [set.seed()]
#' beforehand for reproducibility, the standard convention for
#' [stats::runif()]/[stats::rnorm()] and the rest of R.
#'
#' A typo in `priors`' names raises an error immediately, for consistency
#' with [kupdate()]/[freeze()] elsewhere in this package.
#'
#' A constraint violated by a sampled value (e.g. a negative draw for a
#' strictly positive `length_scale`) raises an error immediately, from
#' inside the relevant [parametrisation]'s `wrap()` -- the same single point
#' of validation used everywhere else in this package (no separate check
#' needed here).
#'
#' @param kernel A kernel object.
#' @param priors A named list, each element a length-2 numeric vector:
#'   `c(min, max)` for `distribution = stats::runif` (the default), or
#'   `c(mean, sd)` for `distribution = stats::rnorm`.
#' @param distribution A sampling function with signature `(n, p1, p2)`,
#'   e.g. [stats::runif()] or [stats::rnorm()]. Defaults to [stats::runif()].
#'
#' @return A new kernel object with the requested hyperparameters resampled.
#' @export
#' @examples
#' set.seed(1)
#' k <- se_kernel(length_scale = 1)
#'
#' k2 <- sample_hps(k, priors = list(length_scale = c(0.1, 5)))
#' get_params(k2)
#'
#' # Normal priors (mean, sd) instead of uniform (min, max):
#' k3 <- sample_hps(k, priors = list(length_scale = c(2, 0.5)), distribution = stats::rnorm)
#'
#' # Frozen hyperparameters are skipped:
#' k4 <- freeze(k, "length_scale")
#' k5 <- sample_hps(k4, priors = list(length_scale = c(0.1, 5)))
#' identical(k4$length_scale, k5$length_scale)
sample_hps <- function(kernel, priors, distribution = stats::runif) {
  if (length(priors) == 0) {
    return(kernel)
  }

  available <- names(flatten_params(kernel))
  missing <- setdiff(names(priors), available)
  if (length(missing) > 0) {
    stop(
      "No such hyperparameter(s) in this kernel: ", paste(missing, collapse = ", "),
      ". Available: ", paste(unique(available), collapse = ", "), ".",
      call. = FALSE
    )
  }

  .sample_hps_recurse(kernel, priors, distribution)
}

.sample_hps_recurse <- function(kernel, priors, distribution) {
  for (name in names(kernel)) {
    field <- kernel[[name]]
    if (inherits(field, "wrapped_param") && name %in% names(priors) && !isTRUE(field$frozen)) {
      p <- priors[[name]]
      sampled <- distribution(length(field$wrapped), p[1], p[2])
      kernel[[name]] <- wrapped_param(sampled, field$parametrisation, field$frozen)
    }
  }
  if (!is.null(kernel$subkernels)) {
    kernel$subkernels <- lapply(kernel$subkernels, .sample_hps_recurse,
      priors = priors, distribution = distribution
    )
  }
  kernel
}
