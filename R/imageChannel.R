#' Image one channel
#'
#' @param x A targeted intensity dataframe (rows: pixels, cols: features)
#' @param coords A spatial coordinates dataframe.
#' @param channel_number The index of the channel to be imaged.
#' @param quantile_lim A parameter that thresholds the maximum intensity values. Default value is 99%. This means that all pixel intensities greater than the 99th percentile are reduced to that of the 99th percentile.
#'
#' @return An image.
#'
#' @import imager
#' @export
#'
#' @examples imageChannel(IntensityDF, coords, channel_number = 4)

imageChannel <- function(final, channel_number = 1, quantile_lim = 0.99) {

  x <- final$IntensityDF
  coords <- final$SpatialCoords

  matrix_image <- matrix(NA, ncol = max(coords$y), nrow = max(coords$x))
  df <- x

  # intensity vector
  x <- x[, channel_number]
  # get intensity value for 99th percentile most intense pixels
  x_max <- quantile(x, probs = quantile_lim)
  # set all pixel values greater than x_max to that of x_max
  x[x > x_max] <- x_max
  matrix_image[as.matrix(coords)] <- x
  # convert matrix to image
  matrix_image <- as.cimg(matrix_image - min(matrix_image, na.rm = TRUE))

  # add colour channels to the image
  matrix_image <- add.color(matrix_image, simple = TRUE)

  # set red and blue channels to zero to get only green
  R(matrix_image) <- 0
  B(matrix_image) <- 0

  plot(matrix_image, main = colnames(df)[channel_number])

}
