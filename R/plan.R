# Layer between evaluate() and evaluate.batch_kernel() (R/wrappers.R): for K
# kernels sharing the same tree shape evaluated on the same X1/X2 (the
# batch_kernel case when batch_over_inputs = FALSE), every stationary/
# dot-product leaf's distance matrix is identical across the K slices --
# only shape() depends on the (per-slice) hyperparameter values. precompute()
# builds a plan once (a tree parallel to the kernel's, with each cacheable
# leaf's distance already computed); eval_plan() re-walks it for a new set
# of hyperparameter values, paying only the shape()/arithmetic cost.
#
# Ported from the validated prototype dev/ameliorations/02-cached-distance.R
# + 03-shape-broadcast.R (see dev/ameliorations/README.md, section C) --
# measured x2 to x6 depending on scenario, and end-to-end in
# dev/demo/magma/magma-demo-3 (Magma/MagmaClust training loops).
#
# Any node the plan doesn't know how to cache (block_kernel,
# input_specific_param_kernel, affine_kernel -- whose distance depends on
# the `offset` hyperparameter, ...) becomes a "fallback" leaf: eval_plan()
# just calls evaluate() on that subtree directly, exactly as today. This is
# per-subtree, not all-or-nothing -- a sum of one cacheable and one
# uncacheable kernel still accelerates the cacheable half.

#' Precompute a distance-based evaluation plan for a kernel tree
#'
#' @description
#' Builds a plan for `kernel`: a tree parallel to its structure where every
#' leaf whose distance matrix does not depend on its hyperparameters (any
#' stationary or dot-product kernel with a [shape()] method, except
#' [affine_kernel()] whose distance depends on `offset`) has that distance
#' precomputed once. [eval_plan()] then re-evaluates the covariance for any
#' kernel sharing this exact tree structure (typically the same kernel after
#' [kupdate()], which never changes structure) by applying only [shape()] to
#' the cached distance -- skipping the distance recomputation entirely.
#'
#' Any node the plan cannot cache (`block_kernel`,
#' `input_specific_param_kernel`, `affine_kernel`, a custom kernel without a
#' `shape()` method, ...) is recorded as a fallback leaf: [eval_plan()] calls
#' [evaluate()] on that subtree directly, exactly as it would without a
#' plan. Missing/infinite values in `x1`/`x2` are detected once here (same
#' point-level rule as [evaluate()]) and applied by [eval_plan()] per leaf,
#' reproducing [evaluate()]'s existing `na_action` semantics exactly.
#'
#' @param kernel A kernel object.
#' @param x1,x2 Inputs, as for [evaluate()].
#' @return An object of class `kernel_plan`.
#' @seealso [eval_plan()]
#' @export
#' @examples
#' k <- se_kernel(length_scale = 1)
#' plan <- precompute(k, c(0, 1, 2))
#' eval_plan(plan, se_kernel(length_scale = 2))
precompute <- function(kernel, x1, x2 = NULL) {
  # Captured before x1/x2 are coerced/transformed, same rule evaluate()
  # itself uses (see ?evaluate) -- eval_plan() defaults to this when its
  # caller doesn't override `.self_covariance` explicitly, so that calling
  # eval_plan() directly (without going through evaluate.batch_kernel(),
  # which always resolves and forwards it) still masks NA correctly.
  is_self_covariance <- is.null(x2) || identical(x1, x2)

  X1_raw <- as_matrix_input(x1)
  X2_raw <- if (is.null(x2)) X1_raw else as_matrix_input(x2)

  na1 <- apply(X1_raw, 1, function(row) any(!is.finite(row)))
  na2 <- apply(X2_raw, 1, function(row) any(!is.finite(row)))
  X1_clean <- X1_raw
  X2_clean <- X2_raw
  X1_clean[!is.finite(X1_clean)] <- 0
  X2_clean[!is.finite(X2_clean)] <- 0

  plan <- .build_plan(kernel, X1_clean, X2_clean, X1_raw, X2_raw)
  structure(
    list(plan = plan, n1 = nrow(X1_raw), n2 = nrow(X2_raw), na1 = na1, na2 = na2,
         is_self_covariance = is_self_covariance),
    class = "kernel_plan"
  )
}

#' @keywords internal
.build_plan <- function(kernel, X1, X2, X1_raw, X2_raw) {
  cls <- class(kernel)[1]
  if (inherits(kernel, "sum_kernel") || inherits(kernel, "product_kernel")) {
    return(list(
      node = "operator",
      class = cls,
      op = if (inherits(kernel, "sum_kernel")) "+" else "*",
      children = lapply(kernel$subkernels, .build_plan,
                         X1 = X1, X2 = X2, X1_raw = X1_raw, X2_raw = X2_raw)
    ))
  }
  if ((inherits(kernel, "stationary_kernel") || inherits(kernel, "dot_product_kernel")) &&
        has_s3_method("shape", kernel) &&
        !inherits(kernel, "affine_kernel")) {
    if (inherits(kernel$engine, "diagonal_engine")) {
      if (nrow(X1) != nrow(X2)) {
        # Let evaluate()/engine_eval.diagonal_engine() raise its own clear
        # error message rather than duplicating that check here.
        return(list(node = "fallback", class = cls, X1 = X1_raw, X2 = X2_raw))
      }
      return(list(node = "shape", class = cls, d = diag_distance(kernel$distance_function, X1, X2)))
    }
    return(list(node = "shape", class = cls, d = distance_matrix(kernel$distance_function, X1, X2)))
  }
  if (inherits(kernel, "variance_kernel") || inherits(kernel, "constant_kernel")) {
    return(list(node = "flat", class = cls, nrow = nrow(X1), ncol = nrow(X2)))
  }
  list(node = "fallback", class = cls, X1 = X1_raw, X2 = X2_raw)
}

#' Re-evaluate a kernel against a precomputed plan
#'
#' @description
#' Re-walks `plan` (built by [precompute()]), applying [shape()] to each
#' cached distance instead of recomputing it -- `kernel` must share the
#' exact tree structure of the kernel the plan was built from (same classes,
#' same nesting; only hyperparameter values may differ, e.g. after
#' [kupdate()]). `na_action`/`mask_variance`/`.self_covariance` reproduce
#' [evaluate()]'s missing-value handling exactly, applied independently per
#' leaf (matching how a sum/product of kernels already combines
#' independently NA-handled sub-results today).
#'
#' @param plan A `kernel_plan`, from [precompute()].
#' @param kernel A kernel object with the same tree structure as the one
#'   `plan` was built from.
#' @param na_action One of `"propagate"` or `"mask"`, as for [evaluate()].
#' @param mask_variance As for [evaluate()].
#' @param ... Forwarded to [evaluate()] for any fallback (uncached) subtree.
#' @param .self_covariance As for [evaluate()]. Defaults to whether the
#'   `x1`/`x2` originally passed to [precompute()] were themselves a
#'   self-covariance call (`x2` omitted or identical to `x1`) -- the same
#'   rule [evaluate()] applies, so calling `eval_plan()` directly (without
#'   going through [batch_kernel()], which always resolves and forwards
#'   this explicitly) still masks `NA` correctly under `na_action = "mask"`.
#' @return The covariance matrix (or vector, for a `diagonal_engine` leaf),
#'   identical to `evaluate(kernel, x1, x2)` for the `x1`/`x2` `plan` was
#'   built from.
#' @seealso [precompute()]
#' @export
eval_plan <- function(plan, kernel, na_action = "propagate", mask_variance = 1, ...,
                       .self_covariance = NULL) {
  is_self_covariance <- if (is.null(.self_covariance)) {
    isTRUE(plan$is_self_covariance)
  } else {
    isTRUE(.self_covariance)
  }
  .eval_plan_node(plan$plan, kernel, plan$na1, plan$na2, na_action, mask_variance,
                   is_self_covariance, ...)
}

#' @keywords internal
.eval_plan_node <- function(node, kernel, na1, na2, na_action, mask_variance,
                             is_self_covariance, ...) {
  cls <- class(kernel)[1]
  if (!identical(cls, node$class)) {
    stop(
      "eval_plan(): `kernel` does not match the tree structure `plan` was built from ",
      "(expected a `", node$class, "` node here, got `", cls, "`). `kernel` must share ",
      "the exact structure (same classes, same nesting) as the kernel passed to ",
      "precompute() -- only hyperparameter values may differ (e.g. after kupdate()).",
      call. = FALSE
    )
  }
  switch(
    node$node,
    operator = {
      left <- .eval_plan_node(node$children[[1]], kernel$subkernels[[1]], na1, na2,
                               na_action, mask_variance, is_self_covariance, ...)
      right <- .eval_plan_node(node$children[[2]], kernel$subkernels[[2]], na1, na2,
                                na_action, mask_variance, is_self_covariance, ...)
      .check_combinable_shapes(left, right, node$class)
      if (node$op == "+") left + right else left * right
    },
    shape = .apply_leaf_na(shape(kernel, node$d), na1, na2, na_action, mask_variance,
                            is_self_covariance),
    flat = {
      value <- if (inherits(kernel, "variance_kernel")) {
        unwrap_param(kernel$variance)
      } else {
        unwrap_param(kernel$value)
      }
      .apply_leaf_na(matrix(value, nrow = node$nrow, ncol = node$ncol), na1, na2,
                      na_action, mask_variance, is_self_covariance)
    },
    fallback = evaluate(kernel, node$X1, node$X2, na_action = na_action,
                         mask_variance = mask_variance, ...,
                         .self_covariance = is_self_covariance)
  )
}

#' @keywords internal
.apply_leaf_na <- function(result, na1, na2, na_action, mask_variance, is_self_covariance) {
  if (na_action == "mask") {
    apply_na_mask(result, na1, na2, is_self_covariance, mask_variance)
  } else {
    apply_na_propagate(result, na1, na2)
  }
}
