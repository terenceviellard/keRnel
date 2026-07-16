# Ergonomic helper, made possible by `+`/`*` being uniformly available on
# every kernel (R/operators.R).
#
# dim_additive_kernel() (the other helper from that design section) is
# deferred: it is defined in terms of active_dims_kernel(), not yet
# implemented (see the build plan, step 7).

#' Build the full interaction expansion of a list of kernels
#'
#' @description
#' For kernels `k1, k2, k3`, the full expansion (`max_order` = 3, the
#' default) is:
#'
#' \deqn{k_1 + k_2 + k_3 + k_1 k_2 + k_1 k_3 + k_2 k_3 + k_1 k_2 k_3}
#'
#' i.e. the sum of every product of 1 to `max_order` of the input kernels,
#' taken without repetition -- `2^n - 1` terms for `n` kernels at the
#' default `max_order = length(kernels)`. Set `max_order` lower (e.g. `2`)
#' to limit the expansion to low-order interactions only.
#'
#' @param kernels A list of kernel objects.
#' @param max_order An integer between `1` and `length(kernels)`. Defaults
#'   to `length(kernels)` (the full expansion).
#'
#' @return A single kernel object (a composition of [sum_kernel()] and
#'   [product_kernel()]).
#' @export
#' @examples
#' k1 <- se_kernel(length_scale = 1)
#' k2 <- linear_kernel(slope_var = 1)
#' k3 <- periodic_kernel(length_scale = 1, period = 2)
#' k <- kernel_interactions(list(k1, k2, k3))
#' evaluate(k, c(0, 1, 2))
#'
#' # Limit to pairwise interactions (k1+k2+k3+k1*k2+k1*k3+k2*k3, no triple product)
#' kernel_interactions(list(k1, k2, k3), max_order = 2)
kernel_interactions <- function(kernels, max_order = length(kernels)) {
  n <- length(kernels)
  if (n < 1) {
    stop("`kernels` must contain at least one kernel.", call. = FALSE)
  }
  if (max_order < 1 || max_order > n) {
    stop("`max_order` must be between 1 and length(kernels).", call. = FALSE)
  }

  terms <- unlist(lapply(seq_len(max_order), function(order) {
    lapply(
      utils::combn(seq_len(n), order, simplify = FALSE),
      function(idx) Reduce(`*`, kernels[idx])
    )
  }), recursive = FALSE)

  Reduce(`+`, terms)
}
