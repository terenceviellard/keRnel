#' Parametrisations: transforming a hyperparameter to/from an unconstrained space
#'
#' @description
#' A parametrisation defines a bijective transformation between a
#' hyperparameter's natural (constrained) space and an unconstrained space
#' suitable for optimisation with a generic, unconstrained optimiser (e.g.
#' `optim(method = "BFGS")`).
#'
#' Every kernel hyperparameter is stored internally in its *wrapped*
#' (unconstrained) form, alongside the parametrisation object that knows how
#' to `wrap()` and `unwrap()` it (see [wrapped_param()]).
#'
#' Numerical validation of the constraint (e.g. strict positivity) happens
#' *inside* `wrap()`, not duplicated in every kernel constructor: every code
#' path that sets a wrapped value -- a kernel constructor, `kupdate()` --
#' goes through `wrap()`, so the constraint only needs to be written once.
#'
#' `unwrap_grad()` is the derivative of `unwrap()` with respect to its own
#' argument (the wrapped value) -- the piece [kernel_grad_free()] needs to
#' chain a [kernel_grad()] gradient (natural space) into the free/wrapped
#' space [get_free_params()]/[set_free_params()] work in:
#' `d(unwrap(p, w))/dw`.
#'
#' @param parametrisation A parametrisation object, e.g. one created by
#'   [log_exp_parametrisation()].
#' @param x A numeric value or vector: in the constrained (natural) space
#'   when calling `wrap()`, in the unconstrained (wrapped) space when calling
#'   `unwrap()`/`unwrap_grad()`.
#'
#' @return A numeric value or vector, in the other space (`wrap()`/
#'   `unwrap()`), or the local derivative (`unwrap_grad()`).
#' @name parametrisation
NULL

#' @rdname parametrisation
#' @export
wrap <- function(parametrisation, x) {
  UseMethod("wrap")
}

#' @rdname parametrisation
#' @export
unwrap <- function(parametrisation, x) {
  UseMethod("unwrap")
}

#' @rdname parametrisation
#' @export
unwrap_grad <- function(parametrisation, x) {
  UseMethod("unwrap_grad")
}

#' Reject non-numeric/non-finite input before any parametrisation transform
#'
#' @description
#' Shared first check for every `wrap.*` method, called before any
#' constraint-specific logic (e.g. positivity). Without it, a non-numeric or
#' non-finite value either slips through silently (`Inf <= 0` is `FALSE`, so
#' a positivity check alone lets `Inf` through) or fails deep inside a
#' comparison/transform with a confusing low-level error (`missing value
#' where TRUE/FALSE needed`, `non-numeric argument to mathematical
#' function`) instead of a clear one naming the parametrisation.
#'
#' @param x The value passed to `wrap()`.
#' @param parametrisation_name The parametrisation's class name, for the
#'   error message.
#' @keywords internal
.validate_finite_numeric <- function(x, parametrisation_name) {
  if (!is.numeric(x) || length(x) == 0 || anyNA(x) || any(!is.finite(x))) {
    stop(
      "`", parametrisation_name, "` requires finite numeric value(s), got ",
      if (is.null(x)) "NULL" else paste0("`", class(x)[1], "`"), ".",
      call. = FALSE
    )
  }
}

#' Identity parametrisation
#'
#' @description
#' Leaves the value unchanged in both directions. Use this for
#' hyperparameters that have no constraint at all, e.g. the additive
#' `value` of `constant_kernel()` (which may be negative).
#'
#' @return An object of class `identity_parametrisation`/`parametrisation`.
#' @export
#' @examples
#' p <- identity_parametrisation()
#' wrap(p, -3.5)
#' unwrap(p, -3.5)
identity_parametrisation <- function() {
  structure(list(), class = c("identity_parametrisation", "parametrisation"))
}

#' @export
wrap.identity_parametrisation <- function(parametrisation, x) {
  .validate_finite_numeric(x, "identity_parametrisation")
  x
}

#' @export
unwrap.identity_parametrisation <- function(parametrisation, x) x

#' @export
unwrap_grad.identity_parametrisation <- function(parametrisation, x) {
  0 * x + 1
}

#' Log-exp parametrisation (strictly positive constraint)
#'
#' @description
#' Maps a strictly positive value to the real line via `log()`, and back via
#' `exp()`. This is the default parametrisation for hyperparameters that
#' must remain strictly positive, e.g. `length_scale` or `variance`.
#'
#' @return An object of class `log_exp_parametrisation`/`parametrisation`.
#' @export
#' @examples
#' p <- log_exp_parametrisation()
#' wrap(p, 1)
#' unwrap(p, 0)
log_exp_parametrisation <- function() {
  structure(list(), class = c("log_exp_parametrisation", "parametrisation"))
}

#' @export
wrap.log_exp_parametrisation <- function(parametrisation, x) {
  .validate_finite_numeric(x, "log_exp_parametrisation")
  if (any(x <= 0)) {
    stop("`log_exp_parametrisation` requires strictly positive value(s).", call. = FALSE)
  }
  log(x)
}

#' @export
unwrap.log_exp_parametrisation <- function(parametrisation, x) exp(x)

#' @export
unwrap_grad.log_exp_parametrisation <- function(parametrisation, x) exp(x)

#' Softplus parametrisation (strictly positive constraint, alternative to log-exp)
#'
#' @description
#' An alternative mapping for strictly positive values, with a gentler
#' gradient near zero than [log_exp_parametrisation()]. Maps to the real
#' line via `log(exp(x) - 1)`, and back via `log(1 + exp(x))`.
#'
#' @return An object of class `softplus_parametrisation`/`parametrisation`.
#' @export
#' @examples
#' p <- softplus_parametrisation()
#' unwrap(p, wrap(p, 3.3))
softplus_parametrisation <- function() {
  structure(list(), class = c("softplus_parametrisation", "parametrisation"))
}

#' @export
wrap.softplus_parametrisation <- function(parametrisation, x) {
  .validate_finite_numeric(x, "softplus_parametrisation")
  if (any(x <= 0)) {
    stop("`softplus_parametrisation` requires strictly positive value(s).", call. = FALSE)
  }
  log(exp(x) - 1)
}

#' @export
unwrap.softplus_parametrisation <- function(parametrisation, x) log(1 + exp(x))

#' @export
unwrap_grad.softplus_parametrisation <- function(parametrisation, x) stats::plogis(x)

#' Bounded parametrisation (interval constraint)
#'
#' @description
#' Maps a value strictly within `(lower, upper)` to the real line via a
#' logit transform, and back via a sigmoid (logistic) transform.
#'
#' @param lower,upper Numeric bounds, with `lower < upper`.
#' @return An object of class `bounded_parametrisation`/`parametrisation`.
#' @importFrom stats plogis qlogis
#' @export
#' @examples
#' p <- bounded_parametrisation(lower = 0, upper = 10)
#' unwrap(p, wrap(p, 4))
bounded_parametrisation <- function(lower, upper) {
  if (any(lower >= upper)) {
    stop("`lower` must be strictly less than `upper`.", call. = FALSE)
  }
  structure(
    list(lower = lower, upper = upper),
    class = c("bounded_parametrisation", "parametrisation")
  )
}

#' @export
wrap.bounded_parametrisation <- function(parametrisation, x) {
  .validate_finite_numeric(x, "bounded_parametrisation")
  lower <- parametrisation$lower
  upper <- parametrisation$upper
  if (any(x <= lower) || any(x >= upper)) {
    stop(
      "`bounded_parametrisation` requires value(s) strictly within (",
      lower, ", ", upper, ").",
      call. = FALSE
    )
  }
  stats::qlogis((x - lower) / (upper - lower))
}

#' @export
unwrap.bounded_parametrisation <- function(parametrisation, x) {
  lower <- parametrisation$lower
  upper <- parametrisation$upper
  lower + (upper - lower) * stats::plogis(x)
}

#' @export
unwrap_grad.bounded_parametrisation <- function(parametrisation, x) {
  (parametrisation$upper - parametrisation$lower) * stats::dlogis(x)
}

#' Chain of parametrisations
#'
#' @description
#' Composes several parametrisations into one: `wrap()` applies each
#' parametrisation in order, `unwrap()` applies them in reverse order. Useful
#' to combine constraints, e.g. a positive value additionally bounded above.
#'
#' @param ... Parametrisation objects, applied in the order given to `wrap()`.
#' @return An object of class `parametrisation_chain`/`parametrisation`.
#' @export
#' @examples
#' chain <- parametrisation_chain(log_exp_parametrisation())
#' unwrap(chain, wrap(chain, 2))
parametrisation_chain <- function(...) {
  parametrisations <- list(...)
  structure(
    list(parametrisations = parametrisations),
    class = c("parametrisation_chain", "parametrisation")
  )
}

#' @export
wrap.parametrisation_chain <- function(parametrisation, x) {
  for (p in parametrisation$parametrisations) x <- wrap(p, x)
  x
}

#' @export
unwrap.parametrisation_chain <- function(parametrisation, x) {
  for (p in rev(parametrisation$parametrisations)) x <- unwrap(p, x)
  x
}

#' @export
unwrap_grad.parametrisation_chain <- function(parametrisation, x) {
  grad <- 1
  for (p in rev(parametrisation$parametrisations)) {
    grad <- grad * unwrap_grad(p, x)
    x <- unwrap(p, x)
  }
  grad
}
