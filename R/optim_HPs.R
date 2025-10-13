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
  input <- as.matrix(unique(input))

  # Calcul de la matrice de covariance
  cov <- pairwise_kernel(kern, input, input) + post_cov

  # Inversion stable avec jitter
  inv <- chol_inv_jitter(cov, pen_diag = pen_diag)

  # Gestion de la moyenne
  if (length(mean) == 1) {
    # Si mean est un scalaire, créer un vecteur de la bonne taille
    n <- nrow(input)
    mean_vec <- rep(mean, n)
  } else {
    # Sinon utiliser mean directement (doit avoir la bonne longueur)
    mean_vec <- mean
    if (length(mean_vec) != nrow(input)) {
      stop("La longueur de 'mean' doit correspondre au nombre de lignes de db$Input")
    }
  }

  # Calcul de la vraisemblance
  log_likelihoods <- dmnorm(db$Output, mean_vec, inv, log = TRUE)
  neg_sum_log_likelihood <- -sum(log_likelihoods)

  # Pour les noyaux additifs, ajouter une pénalité de parcimonie (L1) sur les variances
  if (is(kern, "AdditiveKernel")) {
    hps <- gt_HPs(kern)
    variance_params <- hps[grepl("variance", names(hps))]
    l1_penalty <- 0.1 * sum(abs(unlist(variance_params)))  # λ=0.1 pour favoriser la parcimonie
    neg_sum_log_likelihood <- neg_sum_log_likelihood + l1_penalty
  }

  return(neg_sum_log_likelihood)
}

#' @noRd
gr_sum_logGaussian <- function(hp, db, mean, kern, post_cov, pen_diag) {
  kern <- set_hyperparameters(kern, hp)
  list_hp <- get_hyperparameter_names(kern)
  output <- db$Output
  input <- db$Input
  input <- as.matrix((unique(input))

  # Calcul de la matrice de covariance et son inverse
  cov <- pairwise_kernel(kern, input, input) + post_cov
  inv <- chol_inv_jitter(cov, pen_diag = pen_diag)

  # Gestion de la moyenne
  if (length(mean) == 1) {
    n <- nrow(input)
    mean_vec <- rep(mean, n)
  } else {
    mean_vec <- mean
  }

  # Termes communs pour le gradient
  prod_inv <- inv %*% (output - mean_vec)
  common_term <- prod_inv %*% t(prod_inv) - inv

  # Calcul du gradient pour chaque hyperparamètre
  grad <- numeric(length(list_hp))
  names(grad) <- list_hp

  for (i in seq_along(list_hp)) {
    hp_name <- list_hp[i]
    kern_deriv <- kernel_deriv(kern, input, input, hp_name)
    grad[i] <- sum(diag(-0.5 * (common_term %*% kern_deriv)))

    # Pour les noyaux additifs, ajouter la dérivée de la pénalité L1
    if (is(kern, "AdditiveKernel") && grepl("variance", hp_name)) {
      # Extraire la valeur actuelle du paramètre
      current_value <- unlist(gt_HPs(kern))[[hp_name]]
      # Dérivée de |θ| est sign(θ)
      sign_term <- ifelse(current_value > 0, 1, -1) * 0.1  # λ=0.1
      grad[i] <- grad[i] + sign_term
    }
  }

  return(grad)
}



#' Optimize Hyperparameters for a kernel (with additive kernel support)
#'
#' @param hp A vector of initial hyperparameters to be optimized.
#' @param db The dataset used for optimization.
#' @param mean The mean vector (must have length equal to nrow(db$Input))
#' @param kern The kernel function
#' @param post_cov The posterior covariance function.
#' @param pen_diag A penalty term added to the diagonal of the covariance matrix for numerical stability.
#' @param verbose A logical value indicating whether to return the full optimization result or just the optimized parameters.
#' @param max_iter Maximum number of iterations. Defaults to 1000.
#' @param lambda L1 regularization parameter for additive kernels. Defaults to 0.1.
#'
#' @return If `verbose` is FALSE, a vector of optimized hyperparameters; otherwise, the full result from the `optim` function.
#' @export
optim_hp_L_BFGS_B <- function(hp, db, mean, kern, post_cov, pen_diag=1e-6, verbose = FALSE, max_iter = 1000, lambda = 0.1) {
  # Stocker les paramètres dans l'environnement des fonctions
  objective <- function(hp) {
    sum_logGaussian(hp, db, mean, kern, post_cov, pen_diag)
  }

  gradient <- function(hp) {
    gr_sum_logGaussian(hp, db, mean, kern, post_cov, pen_diag)
  }

  # Définir les bornes pour les hyperparamètres
  lower <- rep(1e-6, length(hp))
  names(lower) <- names(hp)
  upper <- rep(1e6, length(hp))
  names(upper) <- names(hp)

  # Appel à l'optimiseur
  result <- stats::optim(
    par = hp,
    fn = objective,
    gr = gradient,
    method = "L-BFGS-B",
    lower = lower,
    upper = upper,
    control = list(factr = 1e7, maxit = max_iter)
  )

  if (!verbose) {
    return(result$par)
  } else {
    return(result)
  }
}


