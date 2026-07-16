# Tree-algebra operation: cleanly remove a sub-kernel from a composite tree
# (see dev/ameliorations/README.md, item 1e). Without this, removing a term
# means either restructuring the kernel by hand or a numeric hack (e.g.
# kupdate(noise = 1e-12) to approximate "no noise", since a strictly
# positive parametrisation forbids exactly 0) -- drop_terms() does the real
# operation: rebuild the tree without the targeted term(s).

#' Remove every matching sum term from a kernel tree
#'
#' @description
#' Walks `kernel` and rebuilds it with every sub-kernel matching `what`
#' removed. The semantics are deliberately strict rather than permissive:
#'
#' - In a sum, removing a term is the natural operation (`k1 + noise` ->
#'   `k1`) -- this is the intended use case.
#' - In a product, "removing" a factor would mean replacing it with the
#'   identity (1) -- almost never what the caller actually wants (removing
#'   a `variance_kernel()` factor doesn't "remove" anything, it rescales
#'   the kernel). This is refused with an explicit error rather than
#'   guessed.
#' - If the match reaches the root itself, or the only remaining branch
#'   of a wrapper, that is an error too: an empty kernel cannot exist.
#'
#' @param kernel A kernel object (any tree).
#' @param what A kernel S3 class name (e.g. `"white_noise_kernel"`), or a
#'   predicate `function(kernel) -> logical`.
#' @return The rebuilt kernel, with the matched term(s) removed.
#' @export
#' @examples
#' k <- se_kernel(length_scale = 1) + white_noise_kernel(noise = 0.1)
#' drop_terms(k, "white_noise_kernel")
drop_terms <- function(kernel, what) {
  matches <- if (is.function(what)) what else function(k) inherits(k, what)
  out <- .drop_terms_rec(kernel, matches)
  if (is.null(out)) {
    stop(
      "drop_terms(): every term matched -- an empty kernel cannot exist. ",
      "(Did you mean to drop a term from a larger sum?)",
      call. = FALSE
    )
  }
  out
}

#' @keywords internal
.drop_terms_rec <- function(kernel, matches) {
  if (matches(kernel)) {
    return(NULL)
  }
  if (inherits(kernel, "sum_kernel")) {
    left <- .drop_terms_rec(kernel$subkernels[[1]], matches)
    right <- .drop_terms_rec(kernel$subkernels[[2]], matches)
    if (is.null(left)) return(right)
    if (is.null(right)) return(left)
    return(sum_kernel(left, right))
  }
  if (inherits(kernel, "product_kernel")) {
    left <- .drop_terms_rec(kernel$subkernels[[1]], matches)
    right <- .drop_terms_rec(kernel$subkernels[[2]], matches)
    if (is.null(left) || is.null(right)) {
      stop(
        "drop_terms(): the matched kernel is a *factor* of a product. ",
        "Dropping a factor silently rescales the kernel -- refuse to guess. ",
        "Restructure the kernel, or drop the whole product term instead.",
        call. = FALSE
      )
    }
    return(product_kernel(left, right))
  }
  if (!is.null(kernel$subkernels)) {
    kernel$subkernels <- lapply(kernel$subkernels, function(k) {
      out <- .drop_terms_rec(k, matches)
      if (is.null(out)) {
        stop(
          "drop_terms(): matched a sub-kernel inside a `", class(kernel)[1],
          "` wrapper -- dropping through this wrapper is ambiguous (there is ",
          "no general rule for what the wrapper should become without its ",
          "inner kernel). Restructure the kernel instead.",
          call. = FALSE
        )
      }
      out
    })
  }
  kernel
}
