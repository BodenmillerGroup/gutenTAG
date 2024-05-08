#' Image one channel
#'
#' @param x A targeted intensity dataframe (rows: pixels, cols: features) OR output of assignMetapeaks If only dataframe/matrix is provided, spatial coordinates must also be provided using coords.
#' @param coords A spatial coordinates dataframe.
#' @param channel_number The index of the channel to be imaged.
#' @param quantile_lim A parameter that thresholds the maximum intensity values. Default value is 99 percent. This means that all pixel intensities greater than the 99th percentile are reduced to that of the 99th percentile.
#' @param interpolate Perform pixel-wise interpolation.
#' @param axes Plot axes when imaging.
#' @param colna If NA values are present, should they be plotted as black background or as transparent?
#'
#' @return An image.
#'
#' @importFrom imager as.cimg
#' @importFrom imager add.colour
#' @importFrom imager R
#' @importFrom imager B
#' @export
#'

imageChannel <- function(x, coords = NA, channel_number = 1, quantile_lim = 0.99, interpolate = FALSE, axes = FALSE, colna = "black") {

  # validity checks
  .valid.imageChannel(x, coords, channel_number, quantile_lim, interpolate, axes, colna)

  # if input is just a matrix or dataframe (eg. if PCA/NMF) do nothing
  if(is.data.frame(x) | is.matrix(x)){
    df <- x
    coords <- coords
  # otherwise treat as list
  }else{
    df <- x$IntensityDF
    coords <- x$SpatialCoords
  }

  # create matrix of NAs
  matrix_image <- matrix(NA, ncol = max(coords$y), nrow = max(coords$x))

  # intensity vector
  x <- df[, channel_number]
  # get intensity value for 99th percentile most intense pixels
  #x_max <- quantile(x, probs = quantile_lim)
  x_max <- quantile(x, probs = quantile_lim, na.rm = TRUE)

  # set all pixel values greater than x_max to that of x_max
  x[x > x_max] <- x_max
  matrix_image[as.matrix(coords)] <- x
  # convert matrix to image
  matrix_image <- imager::as.cimg(matrix_image - min(matrix_image, na.rm = TRUE))

  # add colour channels to the image
  matrix_image <- imager::add.color(matrix_image, simple = TRUE)

  # set red and blue channels to zero to get only green
  R(matrix_image) <- 0
  B(matrix_image) <- 0

  if(colna == "black"){
    bg <- 1
  }

  if(colna == "white"){
    bg <- 0
  }

  plot(matrix_image, main = colnames(df)[channel_number],
       interpolate = interpolate,
       xlim = c(1, max(coords$x)), ylim = c(max(coords$y), 1),
       axes = axes,
       col.na = rgb(0, 0, 0, bg)) # black background for NA


}
