# Opt-in diagnostics, deliberately separate from validate_kernel(). The
# distinction: validate_kernel() runs on every construction, is cheap
# (structural checks only: right fields, right types), and must never be
# skipped. check_kernel_contract() is heavier (evaluates the kernel,
# computes eigenvalues), meant to be run by hand -- typically once, while
# writing a custom kernel (see vignette("b-intermediate"), "Writing your
# own kernel") -- not on every construction.

#' Detect hyperparameter names shared across sibling sub-kernels
#'
#' @description
#' [kupdate()]/[freeze()] intentionally propagate a named update to *every*
#' occurrence of that name anywhere in the kernel tree -- documented, tested
#' behaviour (see [kupdate()]), not a bug. Sharing a hyperparameter name
#' across sub-kernels is often deliberate: [kernel_interactions()], for
#' instance, routinely combines kernels from different families that happen
#' to both use `length_scale`.
#'
#' This function does not error or warn on its own -- it is a diagnostic to
#' call when you *suspect* an accidental collision (e.g. after a typo in a
#' hand-built composite kernel), not a validation gate on every
#' construction. [get_params()]/[get_trainable_params()] already return
#' safely disambiguated names regardless of whether you check this first.
#'
#' @param kernel A kernel object.
#' @return A character vector of hyperparameter names that appear more than
#'   once in `kernel`'s tree (empty if none).
#' @export
#' @examples
#' k <- se_kernel(length_scale = 1) + periodic_kernel(length_scale = 2, period = 3)
#' check_param_names(k) # "length_scale"
#'
#' k2 <- se_kernel(length_scale = 1) + linear_kernel(slope_var = 1)
#' check_param_names(k2) # character(0): no shared name
check_param_names <- function(kernel) {
  nms <- names(flatten_params(kernel))
  unique(nms[duplicated(nms)])
}

#' Run a battery of runtime sanity checks on a kernel
#'
#' @description
#' A heavier, opt-in complement to [validate_kernel()] (which only checks
#' structure, at construction time, and must stay cheap). Meant to be run by
#' hand -- typically once, while writing a custom kernel (see
#' `vignette("b-intermediate")`, "Writing your own kernel") -- to catch
#' mistakes `validate_kernel()` cannot: a formula that isn't symmetric, a
#' covariance matrix that isn't positive semi-definite, an `evaluate()`
#' method that mishandles the self-covariance case, or a hyperparameter name
#' accidentally shared with a sibling sub-kernel.
#'
#' Checks performed, evaluated on a small set of points (`X`, or `n_points`
#' random ones if `X` is omitted):
#' - `symmetric`: `evaluate(kernel, X)` is symmetric.
#' - `positive_semidefinite`: its eigenvalues are `>= -tol` (relative to
#'   their own scale).
#' - `self_covariance_consistent`: `evaluate(kernel, X)` equals
#'   `evaluate(kernel, X, X)` (the self-covariance and explicit-`x2` code
#'   paths must agree).
#' - `kupdate_roundtrip`: reapplying each trainable hyperparameter's current
#'   value via [kupdate()] reproduces the same result -- `NA` if no
#'   hyperparameter has a single, unambiguous current value to round-trip
#'   (i.e. every trainable name is either unique or already shares an equal
#'   value at every occurrence; see `shared_param_names`).
#' - `shared_param_names`: the result of [check_param_names()], for
#'   information -- not itself a pass/fail check (a shared name is often
#'   intentional).
#'
#' A kernel using an engine other than the default dense one (e.g.
#' [diagonal_engine()], whose `evaluate()` returns a vector, not a matrix)
#' skips `symmetric`/`positive_semidefinite` (`NA`): those checks assume a
#' full covariance matrix.
#'
#' @param kernel A kernel object.
#' @param X Points to evaluate at. Defaults to `n_points` random 1-D points.
#' @param n_points Number of random points to use if `X` is omitted.
#'   Defaults to `5`.
#' @param tol Numerical tolerance. Defaults to `1e-6`.
#'
#' @return An object of class `kernel_contract_check` (a list of named
#'   results), printed as a short pass/fail summary.
#' @export
#' @examples
#' check_kernel_contract(se_kernel(length_scale = 1))
check_kernel_contract <- function(kernel, X = NULL, n_points = 5, tol = 1e-6) {
  if (!inherits(kernel, "kernel")) {
    stop("`kernel` must be a `kernel` object.", call. = FALSE)
  }
  X <- if (is.null(X)) matrix(stats::rnorm(n_points), ncol = 1) else as_matrix_input(X)

  K <- evaluate(kernel, X)
  results <- list()

  if (is.matrix(K)) {
    results$symmetric <- isTRUE(all.equal(K, t(K), tolerance = tol, check.attributes = FALSE))
    eig <- eigen(K, symmetric = TRUE, only.values = TRUE)$values
    results$positive_semidefinite <- all(eig >= -tol * max(1, max(abs(eig))))
  } else {
    results$symmetric <- NA
    results$positive_semidefinite <- NA
  }

  K_explicit <- evaluate(kernel, X, X)
  results$self_covariance_consistent <- isTRUE(
    all.equal(K, K_explicit, tolerance = tol, check.attributes = FALSE)
  )

  results$kupdate_roundtrip <- .check_kupdate_roundtrip(kernel, X, K, tol)
  results$shared_param_names <- check_param_names(kernel)

  structure(results, class = "kernel_contract_check")
}

#' @keywords internal
.check_kupdate_roundtrip <- function(kernel, X, K, tol) {
  trainable <- Filter(function(p) !isTRUE(p$frozen), flatten_params(kernel))
  by_name <- split(trainable, names(trainable))

  # A name is only safe to round-trip through kupdate() (which propagates
  # to *every* occurrence) if it is unique, or if every occurrence already
  # shares the same current value -- otherwise "reapplying" one occurrence's
  # value would silently change the others, which is not a round-trip at
  # all, just a different kernel.
  safe_updates <- list()
  for (nm in names(by_name)) {
    vals <- lapply(by_name[[nm]], unwrap_param)
    agree <- length(vals) == 1 ||
      all(vapply(vals[-1], function(v) isTRUE(all.equal(v, vals[[1]])), logical(1)))
    if (agree) {
      safe_updates[[nm]] <- vals[[1]]
    }
  }

  if (length(safe_updates) == 0) {
    return(NA)
  }
  k2 <- do.call(kupdate, c(list(kernel), safe_updates))
  isTRUE(all.equal(K, evaluate(k2, X), tolerance = tol, check.attributes = FALSE))
}

#' @export
print.kernel_contract_check <- function(x, ...) {
  fmt_check <- function(name, value) {
    status <- if (is.na(value)) "  n/a  " else if (isTRUE(value)) "  OK   " else "FAILED "
    paste0("[", status, "] ", name)
  }
  cat("check_kernel_contract():\n")
  cat(fmt_check("symmetric", x$symmetric), "\n")
  cat(fmt_check("positive_semidefinite", x$positive_semidefinite), "\n")
  cat(fmt_check("self_covariance_consistent", x$self_covariance_consistent), "\n")
  cat(fmt_check("kupdate_roundtrip", x$kupdate_roundtrip), "\n")
  if (length(x$shared_param_names) > 0) {
    cat(
      "[ info  ] shared_param_names: ", paste(x$shared_param_names, collapse = ", "),
      " (see ?check_param_names -- often intentional, not necessarily a problem)\n",
      sep = ""
    )
  } else {
    cat("[  OK   ] shared_param_names: none\n")
  }
  invisible(x)
}
