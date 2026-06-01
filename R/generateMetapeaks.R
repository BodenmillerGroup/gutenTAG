#' Metapeak Generation
#'
#' We define a metapeak as a single peak that represents a cluster of peaks that correspond to the same molecular species.
#'
#' Generation of metapeaks can be described in the following manner:
#'
#'
#' \enumerate{
#'   \item Frequencies of each detected peak are counted for each pixel.
#'   \item Gaussian smoothing is performed over these frequencies.
#'   \item Clustering (using Watershed algorithm) is performed to segment peaks.
#'   \item Peaks that are present in <1\% (default value of tunable threshold parameter) of pixels are filtered out as part of watershed segmentation boundary.
#'
#' }
#'
#' @param x A Cardinal MSImagingExperiment object after peak detection is performed.
#' @param threshold A threshold for peak detection.
#' @param hist_smooth_factor An optional factor for smoothing the histogram of peak counts differently (default: 1.0).
#' @param fixed.limits A tolerance parameter for setting the metapeak limits to be a fixed value centered around the metapeak max. The total width of the metapeak will be twice the value of this parameter.
#' @param sparsity A sparsity parameter controlling the permissible distance between metapeaks. Higher sparsity values enforce greater spacing (fewer metapeaks detected); lower values allow more closely-spaced metapeaks.
#'
#' @return A list containing information on the location, limits and width of
#'   metapeaks, with a \code{$params} element carrying all pipeline parameters
#'   propagated from upstream steps.
#'
#' @export generateMetapeaks
#'
#' @examples
#' path <- system.file("extdata/Example_data/Example_data.imzML", package = "gutenTAG")
#' panel_path <- system.file("extdata/ref_list.csv", package = "gutenTAG")
#' panel <- readPanel(path = panel_path)
#' raw <- readMSIData(path)
#' pre <- preProcess(raw)
#' peaks <- peakDetection(pre)
#' generateMetapeaks(peaks, hist_smooth_factor = 1)

generateMetapeaks <- function(x, threshold = 0.01, hist_smooth_factor = 1, fixed.limits = NULL, sparsity = 3) {

  upstream_params <- list()
  if (is.list(x) && !is(x, "MSImagingExperiment")) {
    upstream_params <- if (!is.null(x$params)) x$params else list()
    x <- x$peaks
  }

  .valid.generateMetapeaks(x, threshold, fixed.limits, sparsity)

  # get counts of pixels in which peak is detected
  count_df <- countPeaks(x)

  detection_threshold <- length(Cardinal::pixels(x)) * threshold

  # smooth counts
  smooth_counts <- smoothPeakCounts(count_df, hist_smooth_factor)

  # rescale smoothed counts so the pixel-count detection_threshold remains
  # meaningful (Gaussian smoothing lowers the histogram maximum)
  scale_factor <- max(count_df$count) / max(smooth_counts$count)
  smooth_counts$count <- smooth_counts$count * scale_factor

  # generate seeds for segmentatation.
  seed_mz <- generateSeedMz(count_df = smooth_counts, detection_threshold = detection_threshold, sparsity = sparsity)

  # metapeak segmentation
  metapeaks <- estimateMetapeaks(count_df = count_df, smooth_count_df = smooth_counts, seed_mz = seed_mz, detection_threshold = detection_threshold, fixed.limits = fixed.limits)

  metapeaks[["params"]] <- c(
    upstream_params,
    list(
      threshold          = threshold,
      hist_smooth_factor = hist_smooth_factor,
      fixed.limits       = fixed.limits,
      sparsity           = sparsity
    )
  )

  return(metapeaks)
}


