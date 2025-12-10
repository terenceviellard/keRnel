## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(keRnel)

## ----eval=FALSE---------------------------------------------------------------
# ?kernel

## -----------------------------------------------------------------------------
#Input for visualization 
x = matrix(rnorm(100), ncol = 2)


## -----------------------------------------------------------------------------
# Create an object of the SEKernel class.
se_kernel = new("SEKernel")
pretty_print(se_kernel)

## -----------------------------------------------------------------------------
get_hyperparameter_values(se_kernel)

hp = c(12,1)
se_kernel = set_hyperparameters(se_kernel, hp)
pretty_print(se_kernel)

## ----warning=FALSE------------------------------------------------------------
visualize_kernel(se_kernel, x)

## -----------------------------------------------------------------------------
plot_kernel_vs_distance(se_kernel, x)

## -----------------------------------------------------------------------------
# Create an object of the SEKernel class.
RQ_kernel = new("RationalQuadraticKernel")
pretty_print(RQ_kernel)

## -----------------------------------------------------------------------------
get_hyperparameter_values(RQ_kernel)

hp = c(12,1,5)
RQ_kernel = set_hyperparameters(RQ_kernel, hp)
pretty_print(RQ_kernel)

## ----warning=FALSE------------------------------------------------------------
visualize_kernel(RQ_kernel, x)

## -----------------------------------------------------------------------------
plot_kernel_vs_distance(RQ_kernel, x)

## -----------------------------------------------------------------------------
# Create an object of the SEKernel class.
per_kernel = new("PeriodicKernel")
pretty_print(per_kernel)

## -----------------------------------------------------------------------------
get_hyperparameter_values(per_kernel)

hp = c(12,1,5)
per_kernel = set_hyperparameters(per_kernel, hp)
pretty_print(per_kernel)

## ----warning=FALSE------------------------------------------------------------
visualize_kernel(per_kernel, x)

## -----------------------------------------------------------------------------
plot_kernel_vs_distance(per_kernel, x)

## -----------------------------------------------------------------------------
# Create an object of the SEKernel class.
mat_kernel12 = new("MaternKernel12")
pretty_print(mat_kernel12)

## -----------------------------------------------------------------------------
get_hyperparameter_values(mat_kernel12)

hp = c(12)
mat_kernel12 = set_hyperparameters(mat_kernel12, hp)
pretty_print(mat_kernel12)

## ----warning=FALSE------------------------------------------------------------
visualize_kernel(mat_kernel12, x)

## -----------------------------------------------------------------------------
plot_kernel_vs_distance(mat_kernel12, x)

## -----------------------------------------------------------------------------
# Create an object of the SEKernel class.
mat_kernel32 = new("MaternKernel32")
pretty_print(mat_kernel32)

## -----------------------------------------------------------------------------
get_hyperparameter_values(mat_kernel32)

hp = c(12)
mat_kernel32 = set_hyperparameters(mat_kernel32, hp)
pretty_print(mat_kernel32)

## ----warning=FALSE------------------------------------------------------------
visualize_kernel(mat_kernel32, x)

## -----------------------------------------------------------------------------
plot_kernel_vs_distance(mat_kernel32, x)

## -----------------------------------------------------------------------------
# Create an object of the SEKernel class.
mat_kernel52 = new("MaternKernel52")
pretty_print(mat_kernel52)

## -----------------------------------------------------------------------------
get_hyperparameter_values(mat_kernel52)

hp = c(12)
mat_kernel52 = set_hyperparameters(mat_kernel52, hp)
pretty_print(mat_kernel52)

## ----warning=FALSE------------------------------------------------------------
visualize_kernel(mat_kernel52, x)

## -----------------------------------------------------------------------------
plot_kernel_vs_distance(mat_kernel52, x)

## -----------------------------------------------------------------------------
# Create an object of the SEKernel class.
lin_kernel = new("LinearKernel")
pretty_print(lin_kernel)

## -----------------------------------------------------------------------------
get_hyperparameter_values(lin_kernel)

hp = c(12,1,5)
lin_kernel = set_hyperparameters(lin_kernel, hp)
pretty_print(lin_kernel)

## ----warning=FALSE------------------------------------------------------------
visualize_kernel(lin_kernel, x)

## -----------------------------------------------------------------------------
plot_kernel_vs_distance(lin_kernel, x)

## -----------------------------------------------------------------------------
# Create an object of the SEKernel class.
cst_kernel = new("ConstantKernel")
pretty_print(cst_kernel)

## -----------------------------------------------------------------------------
get_hyperparameter_values(cst_kernel)

hp = c(12)
cst_kernel = set_hyperparameters(cst_kernel, hp)
pretty_print(cst_kernel)

## ----warning=FALSE------------------------------------------------------------
visualize_kernel(cst_kernel, x)

## -----------------------------------------------------------------------------
plot_kernel_vs_distance(cst_kernel, x)

## -----------------------------------------------------------------------------
# Create an object of the SEKernel class.
noise_kernel = new("NoiseKernel")
pretty_print(noise_kernel)

## -----------------------------------------------------------------------------
get_hyperparameter_values(noise_kernel)

hp = c(12)
noise_kernel = set_hyperparameters(noise_kernel, hp)
pretty_print(noise_kernel)

## ----warning=FALSE------------------------------------------------------------
visualize_kernel(noise_kernel, x)

## -----------------------------------------------------------------------------
plot_kernel_vs_distance(noise_kernel, x)

## ----eval=FALSE---------------------------------------------------------------
# # Access the documentation for a specific kernel
# ?SEKernel
# ?RationalQuadraticKernel
# # Or browse all available kernels
# help(package = "keRnel")

## -----------------------------------------------------------------------------
linear_k = new("LinearKernel",sigma2_b = 1, sigma2_v = 1, c = 0)
noise_k = new("NoiseKernel",value_c = 0.1)

combined_k = linear_k + noise_k

pretty_print(combined_k)

## ----warning=FALSE------------------------------------------------------------
visualize_kernel(combined_k, x)

## -----------------------------------------------------------------------------
plot_kernel_vs_distance(combined_k, x)

## -----------------------------------------------------------------------------
se_space_k = new("SEKernel",variance_se = 1, length_scale_se = 1)  # Space
periodic_time_k = new("PeriodicKernel",variance_per = 1, length_scale_per = 0.5, period = 1)  # Time

combined_k  = se_space_k * periodic_time_k

pretty_print(combined_k)

## ----warning=FALSE------------------------------------------------------------
visualize_kernel(combined_k, x)

## -----------------------------------------------------------------------------
plot_kernel_vs_distance(combined_k, x)

## -----------------------------------------------------------------------------
# Define a custom kernel class
setClass(
  "MinkowskiKernel",
  contains = "AbstractKernel",
  slots = c(
    p = "numeric",    # Minkowski distance shape parameter
    sigma = "numeric" # lengthscale
  )
)

# Define the pairwise kernel method
setMethod(
  "pairwise_kernel", "MinkowskiKernel",
  function(obj, x, y) {
    dx = outer(rowSums(x^2), rowSums(y^2), FUN = "+") - 2 * tcrossprod(x, y)
    exp(- dx^obj@p/ obj@sigma)
    
  }
)

# Define the hyperparameters method
setMethod(
  "gt_HPs", "MinkowskiKernel",
  function(obj) {
    list(sigma = obj@sigma, p = obj@p)
  }
)

# Define the pretty print method
setMethod(
  "pretty_print", "MinkowskiKernel",
  function(obj) {
    sprintf("MinkowskiKernel(p=%.2f, sigma=%.2f)", obj@p, obj@sigma)
  }
)

# Créer une instance du noyau
minkowski_kernel = new("MinkowskiKernel", p = 2, sigma = 1.0)
pretty_print(minkowski_kernel)

# Compute the kernel matrix
K_custom = pairwise_kernel(minkowski_kernel, x, x)

print(K_custom[12,3])

## ----warning=FALSE------------------------------------------------------------
visualize_kernel(minkowski_kernel, x)

## -----------------------------------------------------------------------------
plot_kernel_vs_distance(minkowski_kernel, x)

