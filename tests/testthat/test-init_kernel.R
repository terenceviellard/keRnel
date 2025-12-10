

# Fonction pour tester ensure_abstract_kernel
test_ensure_abstract_kernel <- function() {
  # Test avec un objet AbstractKernel
  abstract_kernel <- new("AbstractKernel")
  result <- ensure_abstract_kernel(abstract_kernel)
  expect_true(is(result, "AbstractKernel"))

  # Test avec un objet numérique qui doit être enveloppé dans un ConstantKernel
  numeric_value <- 2.5
  result <- ensure_abstract_kernel(numeric_value)
  expect_true(is(result, "ConstantKernel"))
  expect_equal(result@value_c, numeric_value)

  # Test avec un objet non valide
  expect_error(ensure_abstract_kernel("invalid"), "The provided kernel is not an AbstractKernel and cannot be converted because it is not numeric.")
}

# Fonction pour tester get_hyperparameter_names
test_get_hyperparameter_names <- function() {
  # Créer un noyau SEKernel pour tester
  se_kernel <- new("SEKernel", variance_se = 1, length_scale_se = 1)
  hyperparam_names <- get_hyperparameter_names(se_kernel)

  # Vérifions que les noms des hyperparamètres sont corrects
  expect_equal(hyperparam_names, c("variance_se", "length_scale_se"))
}

# Fonction pour tester get_hyperparameter_values
test_get_hyperparameter_values <- function() {
  # Créer un noyau SEKernel pour tester
  se_kernel <- new("SEKernel", variance_se = 2, length_scale_se = 0.5)
  hyperparam_values <- get_hyperparameter_values(se_kernel)

  # Vérifions que les valeurs des hyperparamètres sont correctes
  expect_equal(hyperparam_values, c(variance_se = 2, length_scale_se = 0.5))
}

# Fonction pour tester set_hyperparameters
test_set_hyperparameters <- function() {
  # Créer un noyau SEKernel pour tester
  se_kernel <- new("SEKernel", variance_se = 1, length_scale_se = 1)

  # Mettre à jour les hyperparamètres
  updated_kernel <- set_hyperparameters(se_kernel, c(variance_se = 2, length_scale_se = 0.5))

  # Vérifions que les hyperparamètres ont été mis à jour correctement
  expect_equal(get_hyperparameter_values(updated_kernel), c(variance_se = 2, length_scale_se = 0.5))

  # Test avec un ensemble incomplet de valeurs (non nommées)
  expect_error(set_hyperparameters(se_kernel, c(2, 0.5, 3)), "The number of values provided does not match the number of slots in the kernel.")

  # Test avec un nom d'hyperparamètre non valide
  expect_error(set_hyperparameters(se_kernel, c(variance_se = 2, invalid_param = 0.5)), "The slot invalid_param does not exist in the kernel.")
}

# Fonction pour tester SumKernel
test_SumKernel <- function() {
  # Créer quelques noyaux simples pour le test
  kernel1 <- new("ConstantKernel", value_c = 1)
  kernel2 <- new("ConstantKernel", value_c = 2)

  # Créer un SumKernel
  sum_kernel <- kernel1 +kernel2

  # Tester la méthode pairwise_kernel
  x <- matrix(c(1, 2, 3, 4), nrow=2)
  y <- matrix(c(5, 6, 7, 8), nrow=2)
  result <- pairwise_kernel(sum_kernel, x, y)
  expected_result <- matrix(rep(3, 4), nrow=2)  # 1 + 2 pour chaque élément
  expect_equal(result, expected_result)

  # Tester gt_HPs
  hyperparam_values <- gt_HPs(sum_kernel)
  expect_equal(hyperparam_values, list(value_c =1,value_c= 2))
}

# Fonction pour tester ProductKernel
test_ProductKernel <- function() {
  # Créer quelques noyaux simples pour le test
  kernel1 <- new("ConstantKernel", value_c = 2)
  kernel2 <- new("ConstantKernel", value_c = 3)

  # Créer un ProductKernel
  product_kernel <- new("ProductKernel", kernels = list(kernel1, kernel2))

  # Tester la méthode pairwise_kernel
  x <- matrix(c(1, 2, 3, 4), nrow=2)
  y <- matrix(c(5, 6, 7, 8), nrow=2)
  result <- pairwise_kernel(product_kernel, x, y)
  expected_result <- matrix(rep(6, 4), nrow=2)  # 2 * 3 pour chaque élément
  expect_equal(result, expected_result)

  # Tester gt_HPs
  hyperparam_values <- gt_HPs(product_kernel)
  expect_equal(hyperparam_values, list(value_c = 2,value_c =  3))
}

# Fonction pour tester ConstantKernel
test_ConstantKernel <- function() {
  # Créer un noyau constant avec une valeur spécifique
  constant_kernel <- new("ConstantKernel", value_c = 5)

  # Tester la méthode pairwise_kernel
  x <- matrix(c(1, 2, 3, 4), nrow=2)
  y <- matrix(c(5, 6, 7, 8), nrow=2)
  result <- pairwise_kernel(constant_kernel, x, y)
  expected_result <- matrix(rep(5, 4), nrow=2)
  expect_equal(result, expected_result)

  # Tester gt_HPs
  hyperparam_values <- gt_HPs(constant_kernel)
  expect_equal(hyperparam_values, list(value_c = 5))

  # Tester kernel_deriv
  deriv <- kernel_deriv(constant_kernel, x, y, param = "value_c")
  expected_deriv <- matrix(rep(1, 4), nrow=2)
  expect_equal(deriv, expected_deriv)
}

# Fonction pour tester SEKernel
test_SEKernel <- function() {
  # Créer un noyau SEKernel avec des valeurs spécifiques
  se_kernel <- new("SEKernel", variance_se = 2, length_scale_se = 1)

  # Tester la méthode pairwise_kernel
  x <- matrix(c(0, 1, 2, 3), nrow=2)
  y <- matrix(c(0, 1, 2, 3), nrow=2)
  result <- pairwise_kernel(se_kernel, x, y)
  # Calculer le résultat attendu manuellement pour vérifier
  # Pour ce cas simple, on sait que la diagonale devrait être variance_se car distance est nulle,
  # mais pour une vérification exacte, il faudrait calculer la matrice complète basée sur les distances.
  # Pour un test simple, on peut s'assurer que les valeurs sont non-négatives et raisonnables.
  expect_true(all(result >= 0))

  # Tester gt_HPs
  hyperparam_values <- gt_HPs(se_kernel)
  expect_equal(hyperparam_values, list(variance_se = 2, length_scale_se = 1))

  # Tester kernel_deriv (par rapport à variance_se)
  deriv_variance <- kernel_deriv(se_kernel, x, y, param = "variance_se")
  expect_true(all(deriv_variance >= 0))  # Les dérivées devraient être positives

  # Tester kernel_deriv (par rapport à length_scale_se)
  deriv_length <- kernel_deriv(se_kernel, x, y, param = "length_scale_se")
  # La dérivée par rapport à la longueur d'échelle peut être positive ou négative selon les distances
  expect_true(all(!is.nan(deriv_length)))
}

# Exécuter tous les tests
test_that("ensure_abstract_kernel fonctionne correctement", {
  test_ensure_abstract_kernel()
})

test_that("get_hyperparameter_names fonctionne correctement", {
  test_get_hyperparameter_names()
})

test_that("get_hyperparameter_values fonctionne correctement", {
  test_get_hyperparameter_values()
})

test_that("set_hyperparameters fonctionne correctement", {
  test_set_hyperparameters()
})

test_that("SumKernel fonctionne correctement", {
  test_SumKernel()
})

test_that("ProductKernel fonctionne correctement", {
  test_ProductKernel()
})

test_that("ConstantKernel fonctionne correctement", {
  test_ConstantKernel()
})

test_that("SEKernel fonctionne correctement", {
  test_SEKernel()
})
test_NoiseKernel <- function() {
  # Créer un noyau NoiseKernel avec une valeur spécifique
  noise_kernel <- new("NoiseKernel", value_c = 3)

  # Tester la méthode pairwise_kernel
  x <- matrix(c(1, 2, 3, 4), nrow=2)
  result <- pairwise_kernel(noise_kernel, x, x)
  expected_result <- diag(rep(3, 2))
  expect_equal(result, expected_result)

  # Tester gt_HPs
  hyperparam_values <- gt_HPs(noise_kernel)
  expect_equal(hyperparam_values, list(value_c = 3))

  # Tester kernel_deriv
  deriv <- kernel_deriv(noise_kernel, x, x, param = "value_c")
  expected_deriv <- diag(rep(1, 2))
  expect_equal(deriv, expected_deriv)
}

test_that("NoiseKernel fonctionne correctement", {
  test_NoiseKernel()
})

# Tests supplémentaires pour LinearKernel
test_LinearKernel <- function() {
  # Créer un noyau LinearKernel avec des valeurs spécifiques
  linear_kernel <- new("LinearKernel", sigma2_b = 1, sigma2_v = 2, c = 0)

  # Tester la méthode pairwise_kernel
  x <- matrix(c(1, 2, 3, 4), nrow=2)
  y <- matrix(c(5, 6, 7, 8), nrow=2)
  result <- pairwise_kernel(linear_kernel, x, y)

  # Calcul des produits centrés
  x_centered <- x - linear_kernel@c
  y_centered <- y - linear_kernel@c
  expected_result <- linear_kernel@sigma2_b + linear_kernel@sigma2_v * tcrossprod(x_centered, y_centered)
  expect_equal(result, expected_result)

  # Tester gt_HPs
  hyperparam_values <- gt_HPs(linear_kernel)
  expect_equal(hyperparam_values, list(sigma2_b = 1, sigma2_v = 2, c = 0))

  # Tester kernel_deriv par rapport à sigma2_b
  deriv_sigma2_b <- kernel_deriv(linear_kernel, x, y, param = "sigma2_b")
  expected_deriv_sigma2_b <- matrix(1, nrow = nrow(x), ncol = nrow(y))
  expect_equal(deriv_sigma2_b, expected_deriv_sigma2_b)

}

test_that("LinearKernel fonctionne correctement", {
  test_LinearKernel()
})

# Tests supplémentaires pour RationalQuadraticKernel
test_RationalQuadraticKernel <- function() {
  # Créer un noyau RationalQuadraticKernel avec des valeurs spécifiques
  rq_kernel <- new("RationalQuadraticKernel", variance_rq = 2, length_scale_rq = 1, alpha_rq = 1)

  # Tester la méthode pairwise_kernel
  x <- matrix(c(0, 1, 0, 1), nrow=2)
  y <- matrix(c(0, 1, 0, 1), nrow=2)
  result <- pairwise_kernel(rq_kernel, x, y)

  # Vérification simple que toutes les valeurs sont positives
  expect_true(all(result >= 0))

  # Tester gt_HPs
  hyperparam_values <- gt_HPs(rq_kernel)
  expect_equal(hyperparam_values, list(variance_rq = 2, length_scale_rq = 1, alpha_rq = 1))

  # Tester kernel_deriv par rapport à variance_rq
  deriv_variance <- kernel_deriv(rq_kernel, x, y, param = "variance_rq")
  expect_true(all(deriv_variance >= 0))  # Les dérivées devraient être positives

  # Tester kernel_deriv par rapport à length_scale_rq
  deriv_length <- kernel_deriv(rq_kernel, x, y, param = "length_scale_rq")
  expect_true(all(!is.nan(deriv_length)))
}

test_that("RationalQuadraticKernel fonctionne correctement", {
  test_RationalQuadraticKernel()
})

# Tests supplémentaires pour PeriodicKernel
test_PeriodicKernel <- function() {
  # Créer un noyau PeriodicKernel avec des valeurs spécifiques
  periodic_kernel <- new("PeriodicKernel", variance_per = 1, length_scale_per = 1, period = 1)

  # Tester la méthode pairwise_kernel
  x <- matrix(c(0, 0.5, 1, 1.5), nrow=2)
  y <- matrix(c(0, 0.5, 1, 1.5), nrow=2)
  result <- pairwise_kernel(periodic_kernel, x, y)

  # Vérification simple que toutes les valeurs sont positives
  expect_true(all(result >= 0))

  # Tester gt_HPs
  hyperparam_values <- gt_HPs(periodic_kernel)
  expect_equal(hyperparam_values, list(variance_per = 1, length_scale_per = 1, period = 1))

  # Tester kernel_deriv par rapport à variance_per
  deriv_variance <- kernel_deriv(periodic_kernel, x, y, param = "variance_per")
  expect_true(all(!is.nan(deriv_variance)))  # Les dérivées devraient être positives
}

test_that("PeriodicKernel fonctionne correctement", {
  test_PeriodicKernel()
})

# Tests supplémentaires pour MaternKernel12
test_MaternKernel12 <- function() {
  # Créer un noyau MaternKernel12 avec une valeur spécifique
  matern_kernel <- new("MaternKernel12", length_scale_mat = 1)

  # Tester la méthode pairwise_kernel
  x <- matrix(c(0, 1, 0, 1), nrow=2)
  y <- matrix(c(0, 1, 0, 1), nrow=2)
  result <- pairwise_kernel(matern_kernel, x, y)

  # Vérification simple que toutes les valeurs sont positives
  expect_true(all(result >= 0))

  # Tester gt_HPs
  hyperparam_values <- gt_HPs(matern_kernel)
  expect_equal(hyperparam_values, list(length_scale_mat = 1))

  # Tester kernel_deriv par rapport à length_scale_mat
  deriv_length <- kernel_deriv(matern_kernel, x, y, param = "length_scale_mat")
  expect_true(all(!is.nan(deriv_length)))
}

test_that("MaternKernel12 fonctionne correctement", {
  test_MaternKernel12()
})

# Tests supplémentaires pour MaternKernel32
test_MaternKernel32 <- function() {
  # Créer un noyau MaternKernel32 avec une valeur spécifique
  matern_kernel <- new("MaternKernel32", length_scale_mat = 1)

  # Tester la méthode pairwise_kernel
  x <- matrix(c(0, 1, 0, 1), nrow=2)
  y <- matrix(c(0, 1, 0, 1), nrow=2)
  result <- pairwise_kernel(matern_kernel, x, y)

  # Vérification simple que toutes les valeurs sont positives
  expect_true(all(result >= 0))

  # Tester gt_HPs
  hyperparam_values <- gt_HPs(matern_kernel)
  expect_equal(hyperparam_values, list(length_scale_mat = 1))

  # Tester kernel_deriv par rapport à length_scale_mat
  deriv_length <- kernel_deriv(matern_kernel, x, y, param = "length_scale_mat")
  expect_true(all(!is.nan(deriv_length)))
}

test_that("MaternKernel32 fonctionne correctement", {
  test_MaternKernel32()
})

# Tests supplémentaires pour MaternKernel52
test_MaternKernel52 <- function() {
  # Créer un noyau MaternKernel52 avec une valeur spécifique
  matern_kernel <- new("MaternKernel52", length_scale_mat = 1)

  # Tester la méthode pairwise_kernel
  x <- matrix(c(0, 1, 0, 1), nrow=2)
  y <- matrix(c(0, 1, 0, 1), nrow=2)
  result <- pairwise_kernel(matern_kernel, x, y)

  # Vérification simple que toutes les valeurs sont positives
  expect_true(all(result >= 0))

  # Tester gt_HPs
  hyperparam_values <- gt_HPs(matern_kernel)
  expect_equal(hyperparam_values, list(length_scale_mat = 1))

  # Tester kernel_deriv par rapport à length_scale_mat
  deriv_length <- kernel_deriv(matern_kernel, x, y, param = "length_scale_mat")
  expect_true(all(!is.nan(deriv_length)))
}

test_that("MaternKernel52 fonctionne correctement", {
  test_MaternKernel52()
})

# Charger le package testthat pour exécuter les tests
library(testthat)

# Tests supplémentaires pour NoiseKernel
test_NoiseKernel_more <- function() {
  # Créer un noyau NoiseKernel avec une valeur spécifique
  noise_kernel <- new("NoiseKernel", value_c = 2)

  # Cas 1: Matrices d'entrée de différentes tailles
  x1 <- matrix(c(1, 2, 3), nrow=1)
  y1 <- matrix(c(4, 5, 6), nrow=1)
  result1 <- pairwise_kernel(noise_kernel, x1, y1)
  expected_result1 <- diag(2,1)
  expect_equal(result1, expected_result1)

  # Cas 2: Grande matrice
  x2 <- matrix(1:100, nrow=50)
  y2 <- matrix(101:200, nrow=50)
  result2 <- pairwise_kernel(noise_kernel, x2, y2)
  expected_result2 <- diag(rep(2, 50))
  expect_equal(result2, expected_result2)
}

test_that("NoiseKernel gère différentes tailles de matrices", {
  test_NoiseKernel_more()
})

# Tests supplémentaires pour LinearKernel
test_LinearKernel_more <- function() {
  # Créer un noyau LinearKernel avec des valeurs spécifiques
  linear_kernel <- new("LinearKernel", sigma2_b = 1, sigma2_v = 2, c = 0)

  # Cas 1: Matrices de différentes tailles
  x1 <- matrix(c(1, 2), nrow=1)
  y1 <- matrix(c(3, 4), nrow=1)
  result1 <- pairwise_kernel(linear_kernel, x1, y1)
  x_centered <- x1 - linear_kernel@c
  y_centered <- y1 - linear_kernel@c
  expected_result1 <- linear_kernel@sigma2_b + linear_kernel@sigma2_v * tcrossprod(x_centered, y_centered)
  expect_equal(result1, expected_result1)

  # Cas 2: Grande matrice
  x2 <- matrix(1:100, nrow=50)
  y2 <- matrix(101:200, nrow=50)
  result2 <- pairwise_kernel(linear_kernel, x2, y2)
  x2_centered <- x2 - linear_kernel@c
  y2_centered <- y2 - linear_kernel@c
  expected_result2 <- linear_kernel@sigma2_b + linear_kernel@sigma2_v * tcrossprod(x2_centered, y2_centered)
  expect_equal(result2, expected_result2)


  # Tester kernel_deriv pour sigma2_v
  deriv_sigma2_v <- kernel_deriv(linear_kernel, x1, y1, param = "sigma2_v")
  expected_deriv_sigma2_v <- tcrossprod(x_centered, y_centered)
  expect_equal(deriv_sigma2_v, expected_deriv_sigma2_v)
}

test_that("LinearKernel gère différentes tailles de matrices", {
  test_LinearKernel_more()
})

# Tests supplémentaires pour RationalQuadraticKernel
test_RationalQuadraticKernel_more <- function() {
  # Créer un noyau RationalQuadraticKernel avec des valeurs spécifiques
  rq_kernel <- new("RationalQuadraticKernel", variance_rq = 2, length_scale_rq = 1, alpha_rq = 1)

  # Cas 1: Matrices de différentes tailles
  x1 <- matrix(c(0, 1), nrow=1)
  y1 <- matrix(c(0, 1), nrow=1)
  result1 <- pairwise_kernel(rq_kernel, x1, y1)
  expect_true(all(result1 >= 0))

  # Cas 2: Grande matrice
  x2 <- matrix(1:100, nrow=50)
  y2 <- matrix(101:200, nrow=50)
  result2 <- pairwise_kernel(rq_kernel, x2, y2)
  expect_true(all(result2 >= 0))

}

test_that("RationalQuadraticKernel gère différentes tailles de matrices", {
  test_RationalQuadraticKernel_more()
})

# Tests supplémentaires pour PeriodicKernel
test_PeriodicKernel_more <- function() {
  # Créer un noyau PeriodicKernel avec des valeurs spécifiques
  periodic_kernel <- new("PeriodicKernel", variance_per = 1, length_scale_per = 1, period = 1)

  # Cas 1: Matrices de différentes tailles
  x1 <- matrix(c(0, 0.5), nrow=1)
  y1 <- matrix(c(0, 0.5), nrow=1)
  result1 <- pairwise_kernel(periodic_kernel, x1, y1)
  expect_true(all(result1 >= 0))

  # Cas 2: Grande matrice
  x2 <- matrix(1:100, nrow=50)
  y2 <- matrix(101:200, nrow=50)
  result2 <- pairwise_kernel(periodic_kernel, x2, y2)
  expect_true(all(result2 >= 0))

}

test_that("PeriodicKernel gère différentes tailles de matrices", {
  test_PeriodicKernel_more()
})

# Tests supplémentaires pour MaternKernel12
test_MaternKernel12_more <- function() {
  # Créer un noyau MaternKernel12 avec une valeur spécifique
  matern_kernel <- new("MaternKernel12", length_scale_mat = 1)

  # Cas 1: Matrices de différentes tailles
  x1 <- matrix(c(0, 1), nrow=1)
  y1 <- matrix(c(0, 1), nrow=1)
  result1 <- pairwise_kernel(matern_kernel, x1, y1)
  expect_true(all(result1 >= 0))

  # Cas 2: Grande matrice
  x2 <- matrix(1:100, nrow=50)
  y2 <- matrix(101:200, nrow=50)
  result2 <- pairwise_kernel(matern_kernel, x2, y2)
  expect_true(all(result2 >= 0))

}

test_that("MaternKernel12 gère différentes tailles de matrices", {
  test_MaternKernel12_more()
})

# Tests supplémentaires pour MaternKernel32
test_MaternKernel32_more <- function() {
  # Créer un noyau MaternKernel32 avec une valeur spécifique
  matern_kernel <- new("MaternKernel32", length_scale_mat = 1)

  # Cas 1: Matrices de différentes tailles
  x1 <- matrix(c(0, 1), nrow=1)
  y1 <- matrix(c(0, 1), nrow=1)
  result1 <- pairwise_kernel(matern_kernel, x1, y1)
  expect_true(all(result1 >= 0))

  # Cas 2: Grande matrice
  x2 <- matrix(1:100, nrow=50)
  y2 <- matrix(101:200, nrow=50)
  result2 <- pairwise_kernel(matern_kernel, x2, y2)
  expect_true(all(result2 >= 0))

}

test_that("MaternKernel32 gère différentes tailles de matrices", {
  test_MaternKernel32_more()
})

# Tests supplémentaires pour MaternKernel52
test_MaternKernel52_more <- function() {
  # Créer un noyau MaternKernel52 avec une valeur spécifique
  matern_kernel <- new("MaternKernel52", length_scale_mat = 1)

  # Cas 1: Matrices de différentes tailles
  x1 <- matrix(c(0, 1), nrow=1)
  y1 <- matrix(c(0, 1), nrow=1)
  result1 <- pairwise_kernel(matern_kernel, x1, y1)
  expect_true(all(result1 >= 0))

  # Cas 2: Grande matrice
  x2 <- matrix(1:100, nrow=50)
  y2 <- matrix(101:200, nrow=50)
  result2 <- pairwise_kernel(matern_kernel, x2, y2)
  expect_true(all(result2 >= 0))
}

test_that("MaternKernel52 gère différentes tailles de matrices", {
  test_MaternKernel52_more()
})


