#' Image one channel
#'
#' @param x Output from \code{assignMetapeaks}, or a plain intensity data.frame
#'   (rows: pixels, cols: features). If a plain data.frame, \code{coords} must also be provided.
#' @param coords A spatial coordinates data.frame with columns \code{x} and \code{y}.
#'   Required when \code{x} is a data.frame or matrix.
#' @param channel Character name or integer index of the channel to display. Default \code{1}.
#' @param quantile_threshold Numeric in [0, 1]. Pixel intensities above this quantile are
#'   capped to that value before plotting. Default \code{0.99}.
#' @param palette Viridis palette option passed to \code{scale_fill_viridis_c}. One of
#'   \code{"viridis"}, \code{"magma"}, \code{"plasma"}, \code{"inferno"}, \code{"cividis"},
#'   etc. Default \code{"viridis"}.
#' @param na_colour Colour used for NA pixels. Default \code{"black"}.
#' @param ... Additional arguments passed to \code{ggplot2::theme()}.
#'
#' @return A \code{ggplot} object.
#'
#' @examples
#' rdata_path <- system.file("extdata/Example_data/Example_processed.Rdata", package = "gutenTAG")
#' load(rdata_path)
#' imageChannel(x = results$processed, channel = 1)
#'
#' @importFrom ggplot2 geom_raster coord_equal scale_fill_viridis_c
#' @importFrom rlang .data
#' @export
#'

imageChannel <- function(x, coords = NULL, channel = 1,
                         quantile_threshold = 0.99,
                         palette = "viridis",
                         na_colour = "black",
                         ...) {

  # validity checks
  .valid.imageChannel(x, coords, channel, quantile_threshold, palette, na_colour)

  # input handling — accept assignMetapeaks list or plain df/matrix
  if (is.data.frame(x) || is.matrix(x)) {
    df     <- x
    coords <- coords
  } else {
    df     <- x$IntensityDF
    coords <- x$SpatialCoords
  }

  # channel resolution — accept name or index
  if (is.character(channel)) {
    channel_idx <- which(colnames(df) == channel)
    if (length(channel_idx) == 0) stop("Channel '", channel, "' not found in the data.", call. = FALSE)
  } else {
    channel_idx <- channel
  }
  channel_name <- colnames(df)[channel_idx]

  # quantile thresholding
  intensity <- df[, channel_idx]
  cap <- quantile(intensity, probs = quantile_threshold, na.rm = TRUE)
  intensity[intensity > cap] <- cap

  # build plot dataframe
  plot_df <- data.frame(x = coords$x, y = coords$y, intensity = intensity)

  # build ggplot
  p <- ggplot(plot_df, aes(x = x, y = y, fill = intensity)) +
    geom_raster() +
    coord_equal() +
    scale_fill_viridis_c(option = palette, na.value = na_colour) +
    ggtitle(channel_name) +
    theme_void(base_size = 12) +
    theme(...)

  return(p)

}
