#' @noRd
chol_inv_jitter <- function(mat, pen_diag) {
  diag(mat) <- diag(mat) + pen_diag

  tryCatch(
    {
      chol_mat <- chol(mat)
      inv_mat <- chol2inv(chol_mat)
      inv_mat
    },
    error = function(e) {
      chol_inv_jitter(mat, 10 * pen_diag)
    }
  )
}


#' @noRd
dmnorm <- function(x, mu, inv_Sigma, log = FALSE) {
  # Ensure x is a matrix
  if (is.vector(x)) {
    x <- matrix(x, nrow = 1)
  }

  n <- nrow(x)
  p <- ncol(x)

  # Ensure mu is compatible with x
  if (is.vector(mu)) {
    if (length(mu) != p) {
      stop("The length of mu must match the number of columns in x.")
    }
    mu <- matrix(rep(mu, n), ncol = p, byrow = TRUE)
  } else if (is.matrix(mu)) {
    if (ncol(mu) != p || nrow(mu) != n) {
      stop("The dimensions of mu must match the dimensions of x.")
    }
  } else {
    stop("mu must be a vector or a matrix.")
  }

  # Compute z as the difference between x and mu
  z <- x - mu


  # Ensure dimensions are compatible for multiplication
  if (ncol(z) != nrow(inv_Sigma)) {
    if (nrow(z) == nrow(inv_Sigma)) {
      # Transpose z if necessary
      z <- t(z)
    } else {
      stop("The number of columns in z must match the number of rows in inv_Sigma.")
    }
  }

  logdetS <- try(-determinant(inv_Sigma, logarithm = TRUE)$modulus, silent = TRUE)
  attributes(logdetS) <- NULL

  # Perform the multiplication correctly
  ssq <- sum((z %*% inv_Sigma) * z)

  loglik <- as.vector(-(n * log(2 * pi) + logdetS + ssq) / 2)

  if (log) {
    return(loglik)
  } else {
    return(exp(loglik))
  }
}






#' @noRd
sum_logGaussian <- function(hp, db, mean, kern, post_cov, pen_diag) {
  kern <- set_hyperparameters(kern, hp)
  input <- db$Input
  input <- as.matrix(input)
  cov <- pairwise_kernel(kern, input, input) + post_cov
  inv <- chol_inv_jitter(cov, pen_diag = pen_diag)

  if (length(mean) == 1) {
    mean <- rep(mean, nrow(db))
  }

  log_likelihoods <- dmnorm(db$Output, mean, inv, log = TRUE)
  neg_sum_log_likelihood <- -sum(log_likelihoods)
  return(neg_sum_log_likelihood)
}




#' @noRd
gr_sum_logGaussian <- function(hp, db, mean, kern, post_cov, pen_diag) {
  kern <- set_hyperparameters(kern, hp)
  list_hp <- get_hyperparameter_names(kern)

  output <- db$Output
  input <- db$Input
  input <- as.matrix(input)

  cov <- pairwise_kernel(kern, input, input) + post_cov

  inv <- chol_inv_jitter(cov, pen_diag = pen_diag)

  prod_inv <- inv %*% (output - mean)
  common_term <- prod_inv %*% t(prod_inv) - inv

  floop <- function(hp_name) {
    kern_deriv <- kernel_deriv(kern, input, input, hp_name)
    grad_term <- sum(diag((-0.5 * (common_term %*% kern_deriv))))
  }

  sapply(list_hp, floop)
}

#' Optimize Hyperparameters for a kernel
#'
#' This function optimizes hyperparameters using the L-BFGS-B optimization method.
#'
#' @param hp A vector of initial hyperparameters to be optimized.
#' @param db The dataset used for optimization.
#' @param mean The mean function
#' @param kern The kernel function
#' @param post_cov The posterior covariance function.
#' @param pen_diag A penalty term added to the diagonal of the covariance matrix for numerical stability.
#' @param verbose A logical value indicating whether to return the full optimization result or just the optimized parameters. Defaults to FALSE.
#'
#' @return If `verbose` is FALSE, a vector of optimized hyperparameters; otherwise, the full result from the `optim` function.
#' @export
optim_hp_L_BFGS_B <- function(hp, db, mean, kern, post_cov, pen_diag=1e-6, verbose = FALSE) {
  # Call the optimization function
  result <- stats::optim(
    par = hp,
    fn = function(hp) sum_logGaussian(hp, db, mean, kern, post_cov, pen_diag),
    gr = function(hp) gr_sum_logGaussian(hp, db, mean, kern, post_cov, pen_diag),
    method = "L-BFGS-B",
    control = list(factr = 1e7, maxit = 1000)
  )
  if (!verbose) {
    return(result$par)
  } else {
    return(result)
  }
}
