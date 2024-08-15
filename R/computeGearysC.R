#' Compute Geary's C score for each marker.
#'
#' @param x Output from assignMetapeaks
#' @param verbose Boolean flag for verbose output
#' @param update_correspondence Should the output from computing Geary's C score be added to the Correspondence Matrix (bool)
#'
#' @return A vector of Geary's C score for each marker.
#' @importFrom N2R Knn
#' @export computeGearysC
#'

computeGearysC <- function(x, verbose = FALSE, update_correspondence = FALSE){

  # validity checks
  .valid.computeGearysC(x, verbose, update_correspondence)

  # extract intensity data.
  df <- x$IntensityDF
  correspondence <- x$CorrespondenceMatrix
  coords <- x$SpatialCoords

  # create spatial weight matrix (adjacency matrix on spatial coordinates for rooks case, k = 4)
  spatial_weight_matrix <- Knn(as.matrix(coords), k = 4, verbose = verbose, indexType = "L2")
  spatial_weight_matrix <- as(spatial_weight_matrix, "dgCMatrix")

  geary_vector <- c()
  for (i in 1:ncol(df)){

    X <- matrix(df[, i], ncol = 1)
    X_squared <- X^2

    product_temp_1 <- 2 * sum(spatial_weight_matrix %*% X_squared)
    product_temp_2 <- 2 * t(X) %*% spatial_weight_matrix %*% X

    # define params for Geary C computation
    N <- length(X)
    W <- sum(spatial_weight_matrix)
    Var_X <- var(X) * N

    # compute Geary's C score
    gearys_C <- as.numeric((N - 1) * (product_temp_1 - product_temp_2) / (2 * Var_X * W))

    # append each Geary C score to vector
    geary_vector <- c(geary_vector, gearys_C)

  }

  # optional argument to control if Geary's C column is added to Correspondence matrix
  if(update_correspondence == TRUE){

    correspondence$GearysC <- geary_vector
    x$CorrespondenceMatrix <- correspondence

    return(x)

  }else{
    return(geary_vector)
  }



}
