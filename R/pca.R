#' Principal Component Analysis
#'
#' @param x A targeted intensity dataframe (rows: pixels, cols: features)
#' @param comp The number of principal components to compute.
#' @param scree Plot scree plot (bool)
#' @param plot Scatter plot of first two principal components (bool)
#'
#' @return A PCA
#' @export
#'
#' @examples
#' pca(IntensityDF, comp = 5)

pca <- function(x, comp = 5, scree = FALSE, plot = FALSE){

  library(irlba)

  PCA <- prcomp_irlba(x, n = comp, scale. = TRUE)

  if (scree == TRUE){

    library(ggplot2)

    PCA_variance <- data.frame(Variation_Explained = (PCA$sdev)^2 / sum((PCA$sdev)^2) * 100,
                               Component = 1:comp)

    scree_plot <- ggplot(data = PCA_variance, aes(x = as.factor(Component), y = Variation_Explained, fill = as.factor(Component))) +
      geom_bar(stat = "identity") +
      labs(title = "Scree Plot", x = "Component", y = "Variance Explained (%)") +
      scale_y_continuous(limits = c(0, 100)) +
      theme(legend.position = "none")
    print(scree_plot)

  }

  if (plot == TRUE){

    library(ggplot2)

    PCA_data <- data.frame(PCA$x)

    scatter_plot <- ggplot(data = PCA_data, aes(x = PC1, y = PC2)) +
      geom_point() +
      labs(title = "Scatter plot of first two principal components", x = "PC1", y = "PC2")
    print(scatter_plot)
  }

  return(PCA)

}





