#' Compute Geary's C score for each marker.
#'
#' @param object Output from getIntensityDF.
#'
#' @return A vector of Geary's C score for each marker.
#' @import N2R
#' @export
#'
#' @examples
#' constructGearysC(IntensityDF, coords)

constructGearysC <- function(object, update_correspondence = FALSE){

  #library(N2R)

  df <- final$IntensityDF
  correspondence <- final$CorrespondenceMatrix
  coords <- final$SpatialCoords

  spatial_weight_matrix <- N2R::Knn(as.matrix(coords), k = 4, verbose = TRUE, indexType = "L2")
  spatial_weight_matrix <- as(spatial_weight_matrix, "dgCMatrix")

  geary_vector <- c()
  for (i in 1:ncol(df)){

    x <- matrix(df[, i], ncol = 1)
    x_squared <- x^2

    product_temp_1 <- 2 * sum(spatial_weight_matrix %*% x_squared)
    product_temp_2 <- 2 * t(x) %*% spatial_weight_matrix %*% x

    N <- length(x)
    W <- sum(spatial_weight_matrix)
    Var_X <- var(x) * N

    # compure Geary's C score
    gearys_C <- as.numeric((N - 1) * (product_temp_1 - product_temp_2) / (2 * Var_X * W))

    geary_vector <- c(geary_vector, gearys_C)

  }

  if(update_correspondence == TRUE){

    correspondence$GearysC <- geary_vector
    final$CorrespondenceMatrix <- correspondence

  }

  return(final)

}
