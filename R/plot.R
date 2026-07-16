# Visual complement to kernel_tree()/format() (R/print.R): every demo in
# dev/demo builds a covariance heatmap by hand (expand.grid() + geom_tile());
# plot_kernel() is that one-call version (see dev/ameliorations/README.md,
# item 1f).
#
# aes() below refers to `i`/`j`/`k` as bare column names (standard ggplot2
# NSE) rather than the `.data$i` pronoun: `.data` only resolves inside the
# data mask ggplot2 builds internally, and rebinding it via an explicit
# `ggplot2::.data` (needed since ggplot2 is a Suggests-only dependency here,
# not importable) breaks exactly that mechanism at plot-build time. The
# global-variable declaration below is the standard workaround for R CMD
# check's static analysis of aes()'s NSE in this situation (see
# `?utils::globalVariables`).
utils::globalVariables(c("i", "j", "k"))

#' Heatmap of a kernel's covariance matrix
#'
#' @description
#' Evaluates `kernel` over `x` and renders the resulting covariance matrix
#' as a `ggplot2` heatmap -- a quick visual sanity check of a kernel's
#' structure (bandwidth, periodicity, block structure...), independent of
#' any particular downstream model.
#'
#' @param kernel A kernel object.
#' @param x Evaluation points (a numeric vector or matrix, as for
#'   [evaluate()]).
#' @return A `ggplot` object.
#' @export
#' @examples
#' k <- se_kernel(length_scale = 1)
#' plot_kernel(k, seq(-5, 5, length.out = 50))
plot_kernel <- function(kernel, x) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("plot_kernel() needs the 'ggplot2' package.", call. = FALSE)
  }
  K <- evaluate(kernel, x)
  n <- nrow(K)
  df <- expand.grid(i = seq_len(n), j = seq_len(n))
  df$k <- as.vector(K)
  ggplot2::ggplot(df, ggplot2::aes(i, j, fill = k)) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_viridis_c() +
    ggplot2::scale_y_reverse() +
    ggplot2::coord_fixed() +
    ggplot2::labs(title = format(kernel), x = NULL, y = NULL, fill = "k(x1, x2)")
}
