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
#' @importFrom imager R
#' @importFrom imager G
#' @export
#'

imageChannel <- function(x, coords = NA, channel_number = 1, quantile_lim = 0.99, interpolate = FALSE, axes = FALSE, colna = "black", cmap = NA) {

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
  if (is.na(cmap)) {
    # Simply fill the green channel with the intensity values
    cRamp <- grDevices::colorRamp(c(grDevices::rgb(0, 0, 0, 1), grDevices::rgb(0, 1, 0, 1)), interpolate="linear", space="rgb")
  } else if (cmap %in% c("magma", "inferno", "plasma", "viridis", "cividis", "rocket", "mako", "turbo")) {
    # Use viridisLite to generate a color ramp
    cRamp <- grDevices::colorRamp(viridisLite::virids(256, option = cmap), interpolate="linear", space="rgb")
  } else {
    # Assume that cmap is a color ramp
    cRamp <- cmap
  }

  # scale the intensities from [min(matrix_image), percentile(matrix_image, quantile_lim)] to [0, 1] by clipping and scaling
  min_intensity <- min(matrix_image, na.rm = TRUE)
  max_intensity <- quantile(matrix_image, quantile_lim, na.rm = TRUE)
  scaled_intensities <- (pmin(matrix_image, max_intensity) - min_intensity) %/% (max_intensity - min_intensity)

  # apply the specified color map
  matrix_image <- cRamp(scaled_intensities)

  # determine background color
  if (colna == "black") {
    bg <- 1
  } else if (colna == "white") {
    bg <- 0
  } else {
    stop("colna must be either 'black' or 'white'")
  }

  plot(matrix_image, main = colnames(df)[channel_number],
       interpolate = interpolate,
       xlim = c(1, max(coords$x)), ylim = c(max(coords$y), 1),
       axes = axes,
       col.na = rgb(0, 0, 0, bg)) # black background for NA

}
