#' Compute Geary's C score for each marker.
#'
#' Computes Geary's C spatial autocorrelation statistic for each marker channel
#' in the intensity matrix. Values close to 0 indicate strong positive spatial
#' autocorrelation; values close to 2 indicate negative autocorrelation; a value
#' of 1 indicates no spatial structure. Channels with zero variance (e.g.
#' all-zero after failed peak detection) or all-NA pixels return \code{NA}.
#'
#' @param x A list. Output from \code{\link{assignMetapeaks}}, containing at
#'   minimum \code{IntensityDF} (pixels x markers intensity matrix),
#'   \code{SpatialCoords} (pixels x 2 coordinate data.frame), and
#'   \code{CorrespondenceMatrix}.
#' @param verbose logical. If \code{TRUE}, passes verbose output through to the
#'   internal \code{Knn} nearest-neighbour call. Default \code{FALSE}.
#' @param update_correspondence logical. If \code{TRUE}, appends a \code{GearysC}
#'   column to \code{x$CorrespondenceMatrix} and returns the full updated list.
#'   If \code{FALSE} (default), returns a numeric vector of Geary's C scores,
#'   one per marker channel.
#'
#' @return If \code{update_correspondence = FALSE}, a numeric vector of length
#'   equal to the number of marker channels, with \code{NA} for channels that
#'   have zero or missing variance. If \code{update_correspondence = TRUE},
#'   the input list \code{x} with a \code{GearysC} column appended to
#'   \code{x$CorrespondenceMatrix}.
#'
#' @examples
#' path <- system.file("extdata/Example_data.imzML", package = "gutenTAG")
#' panel_path <- system.file("extdata/ref_list.csv", package = "gutenTAG")
#' panel <- readPanel(path = panel_path)
#' raw <- readMSIData(path)
#' pre <- preProcess(raw, cores = 2)
#' peaks <- peakDetection(pre, cores = 2)
#' metapeaks <- generateMetapeaks(peaks)
#' processed <- assignMetapeaks(metapeaks, pre = pre, refList = panel)
#' computeGearysC(processed, verbose = FALSE)
#'
#' @importFrom N2R Knn
#' @importFrom Matrix t
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

  geary_vector <- numeric(ncol(df))
  for (i in seq_len(ncol(df))){

    X <- matrix(df[, i], ncol = 1)
    X_squared <- X^2

    product_temp_1 <- 2 * sum(spatial_weight_matrix %*% X_squared)
    product_temp_2 <- 2 * t(X) %*% spatial_weight_matrix %*% X

    # define params for Geary C computation
    N <- length(X)
    W <- sum(spatial_weight_matrix)
    # na.rm = TRUE so that pixel-level NAs (e.g. tissue edge pixels) are silently
    # dropped rather than propagating NA into Var_X. This is a deliberate choice —
    # revisit if a warning or error on NA-containing channels is preferred.
    Var_X <- var(X, na.rm = TRUE) * N

    # compute Geary's C score — assign NA for zero-variance channels
    if (is.na(Var_X) || Var_X == 0) {
      gearys_C <- NA_real_
    } else {
      gearys_C <- as.numeric((N - 1) * (product_temp_1 - product_temp_2) / (2 * Var_X * W))
    }

    geary_vector[i] <- gearys_C

  }

  # optional argument to control if Geary's C column is added to Correspondence matrix
  if (update_correspondence) {

    correspondence$GearysC <- geary_vector
    x$CorrespondenceMatrix <- correspondence

    return(x)

  }else{
    return(geary_vector)
  }



}
