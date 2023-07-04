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
  
  PCA <- prcomp_irlba(Final_intensity_matrix_targeted, n = comp)
  
  if (scree == TRUE){
    
    PCA_variance <- data.frame(Variation_Explained = (PCA_simple$sdev)^2 / PCA_simple$totalvar * 100,
                               Component = colnames(PCA_simple$x))
    
    ggplot(data = PCA_variance, aes(x = as.factor(Component), y = Variation_Explained, fill = Component)) +
      geom_bar(stat = "identity") +
      labs(title = "Scree Plot", x = "", y = "Variance Explained") +
      scale_y_continuous(limit = c(0, 100))
  }

}

