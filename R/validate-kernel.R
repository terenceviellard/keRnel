# S3 gives no compile-time guarantee that a kernel of a given class has the
# fields its category requires. validate_kernel() is the explicit,
# conventional substitute: a table of per-class contracts, called once at
# the end of every kernel constructor. The real safety net is the companion
# meta-test in tests/testthat/test-validate-kernel.R, which checks that
# every exported constructor actually produces a conformant object.

.kernel_contracts <- list(
  # Runs for *every* kernel (the "kernel" class is always present in
  # class(kernel)): the layer-1-vs-layer-4 contract from
  # dev/ameliorations/design/01-architecture.Rmd, checked generically via S3 method presence
  # rather than listed per concrete class -- a new kernel is covered
  # automatically, with no entry to remember adding here.
  kernel = function(kernel) {
    is_composite <- inherits(kernel, "operator_kernel") || inherits(kernel, "wrapper_kernel")
    if (!is_composite) {
      if (!has_s3_method("pairwise", kernel) && !has_s3_method("evaluate", kernel)) {
        stop(
          "`", class(kernel)[1], "` implements neither `pairwise()` nor its own ",
          "`evaluate()` method -- every base kernel must implement `pairwise()` ",
          "(see vignette(\"b-intermediate\"), \"Writing your own kernel\").",
          call. = FALSE
        )
      }
      if (!inherits(kernel$engine, "engine")) {
        stop(
          "`", class(kernel)[1], "` must have an `engine` field that is an `engine` object.",
          call. = FALSE
        )
      }
    }
  },
  stationary_kernel = function(kernel) {
    if (!is.function(kernel$distance_function)) {
      stop(
        "A `stationary_kernel` must have a `distance_function` field that is a function.",
        call. = FALSE
      )
    }
  },
  dot_product_kernel = function(kernel) {
    if (!is.function(kernel$distance_function)) {
      stop(
        "A `dot_product_kernel` must have a `distance_function` field that is a function.",
        call. = FALSE
      )
    }
  },
  operator_kernel = function(kernel) {
    if (!is.list(kernel$subkernels) || length(kernel$subkernels) != 2) {
      stop("An `operator_kernel` must have exactly 2 `subkernels`.", call. = FALSE)
    }
    if (!all(vapply(kernel$subkernels, inherits, logical(1), what = "kernel"))) {
      stop("Every `subkernel` of an `operator_kernel` must itself be a `kernel`.", call. = FALSE)
    }
    if (!has_s3_method("evaluate", kernel)) {
      stop("An `operator_kernel` must override `evaluate()` directly.", call. = FALSE)
    }
  },
  wrapper_kernel = function(kernel) {
    if (!is.list(kernel$subkernels) || length(kernel$subkernels) != 1) {
      stop("A `wrapper_kernel` must have exactly 1 `subkernel`.", call. = FALSE)
    }
    if (!inherits(kernel$subkernels[[1]], "kernel")) {
      stop("The `subkernel` of a `wrapper_kernel` must itself be a `kernel`.", call. = FALSE)
    }
    if (!has_s3_method("evaluate", kernel)) {
      stop("A `wrapper_kernel` must override `evaluate()` directly.", call. = FALSE)
    }
  }
)

#' Validate that a kernel object satisfies the contract of its classes
#'
#' @description
#' Call this once at the end of every custom kernel constructor (see
#' `vignette("b-intermediate")`, "Writing your own kernel") -- the same
#' convention every built-in kernel constructor follows. It is the
#' explicit substitute, in an S3 world with no compile-time field
#' guarantees, for the structural checks a statically typed class system
#' would enforce automatically.
#'
#' @param kernel A kernel object.
#' @return `kernel`, invisibly, if every contract holds. Errors otherwise.
#' @export
#' @examples
#' pairwise.my_kernel <- function(kernel, x1, x2) kernel$distance_function(x1, x2)
#' k <- structure(
#'   list(distance_function = euclidean_distance, engine = dense_engine()),
#'   class = c("my_kernel", "stationary_kernel", "kernel")
#' )
#' validate_kernel(k) # passes: pairwise() is implemented, engine is valid
validate_kernel <- function(kernel) {
  for (cls in class(kernel)) {
    contract <- .kernel_contracts[[cls]]
    if (!is.null(contract)) {
      contract(kernel)
    }
  }
  invisible(kernel)
}
