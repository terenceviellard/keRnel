#' @title Visualize a Kernel Matrix
#' @description Plots a heatmap of the kernel matrix computed from input data.
#'
#' @param kernel A kernel object for which to compute the pairwise kernel matrix.
#' @param x A matrix or data frame where each row represents an input data point.
#'
#' @return A ggplot object displaying the kernel matrix.
#' @export
#' @importFrom ggplot2 ggplot aes geom_tile scale_fill_gradient labs theme_minimal theme element_blank element_rect
visualize_kernel <- function(kernel, x) {
  # Calculer la matrice du noyau
  kernel_matrix <- pairwise_kernel(kernel, x, x)

  # Convertir la matrice en un data frame pour ggplot2
  df <- as.data.frame(as.table(kernel_matrix))
  names(df) <- c("X1", "X2", "Value")
  df$X2 <- factor(df$X2, levels = rev(unique(df$X2)))

  # Créer un ggplot
  ggplot(df, aes(x = X1, y = X2, fill = Value)) +
    geom_tile() +
    scale_fill_gradient(low = "white", high = "blue") +
    labs(title = "Kernel Matrix Visualization", x = NULL, y = NULL) +
    theme_minimal() +
    theme(
      axis.text.x = element_blank(),
      axis.text.y = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, size = 1)
    )
}

#' @title Plot Kernel Values vs Distance
#' @description Plots kernel values against the distance between input values.
#'
#' @param kernel A kernel object for which to compute the pairwise kernel values.
#' @param input_values A vector of input values for which to compute distances and kernel values.
#'
#' @return A ggplot object showing the kernel values versus distances.
#' @export
#' @importFrom ggplot2 ggplot aes geom_point geom_line labs theme_minimal
plot_kernel_vs_distance <- function(kernel, input_values) {
  # Assurez-vous que les valeurs sont triées pour une meilleure visualisation
  input_values <- sort(input_values)
  names(input_values) <- paste0("Input", 1:length(input_values)) # Noms pour les entrées

  # Calculer la matrice du noyau
  kernel_matrix <- pairwise_kernel(kernel, as.matrix(input_values), as.matrix(input_values))

  # Calculer la matrice des distances absolues
  distance_matrix <- outer(input_values, input_values, FUN = "-")

  # Convertir les matrices en data frames
  df_kernel <- as.data.frame(as.table(kernel_matrix))
  df_distance <- as.data.frame(as.table(distance_matrix))

  # Fusionner les data frames
  df <- data.frame(
    Input1 = df_kernel$Var1,
    Input2 = df_kernel$Var2,
    KernelValue = df_kernel$Freq,
    Distance = df_distance$Freq
  )

  # Tracer le graphique
  ggplot(df, aes(x = Distance, y = KernelValue)) +
    geom_point() +
    geom_line(color = "blue") +
    labs(title = "Kernel vs Distance", x = "Distance", y = "Kernel Value") +
    theme_minimal()
}
