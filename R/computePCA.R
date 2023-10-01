#' Principal Component Analysis
#'
#' @param x A targeted intensity dataframe (rows: pixels, cols: features)
#' @param comp The number of principal components to compute.
#' @param scree Plot scree plot (bool)
#' @param plot Scatter plot of first two principal components (bool)
#'
#' @return A PCA
#'
#' @importFrom irlba prcomp_irlba
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 geom_bar
#' @importFrom ggplot2 labs
#' @importFrom ggplot2 scale_y_continuous
#' @importFrom ggplot2 theme
#' @export computePCA
#'
#' @examples
#' pca(x$IntensityDF, comp = 5)

computePCA <- function(x, comp = 5, seed = 123, scree = FALSE, plot = FALSE){

  # validity checks
  .valid.pca(x, comp, seed, scree, plot)

  # set seed for irlba pca
  set.seed(seed)
  PCA <- prcomp_irlba(x, n = comp, scale. = TRUE)

  if (scree == TRUE){

    PCA_variance <- data.frame(Variation_Explained = (PCA$sdev)^2 / sum((PCA$sdev)^2) * 100,
                               Component = 1:comp)

    scree_plot <- ggplot(data = PCA_variance, ggplot2::aes(x = as.factor(Component), y = Variation_Explained, fill = as.factor(Component))) +
      geom_bar(stat = "identity") +
      labs(title = "Scree Plot", x = "Component", y = "Variance Explained (%)") +
      scale_y_continuous(limits = c(0, 100)) +
      theme(legend.position = "none")

    print(scree_plot)

  }

  if (plot == TRUE){

    PCA_data <- data.frame(PCA$x)

    scatter_plot <- ggplot(data = PCA_data, aes(x = PC1, y = PC2)) +
      geom_point() +
      labs(title = "Scatter plot of first two principal components", x = "PC1", y = "PC2")

    print(scatter_plot)
  }

  return(PCA)

}





