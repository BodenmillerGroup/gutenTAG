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
#'   \item Peaks that are present in <1\% of pixels are filtered out as part of watershed segmentation boundary.
#'
#' }
#'
#' @param x A Cardinal MSImagingExperiment object after peak detection is performed.
#' @param threshold A threshold for peak detection.
#' @param hist_smooth_factor An optional factor for smoothing the histogram of peak counts differently (default: 1.0).
#'
#' @return A list containing information on the location, delimitations and width of metapeaks
#'
#' @import Cardinal
#' @export generateMetapeaks
#'
#' @examples
#' path <- system.file("extdata/Example_data.imzML", package = "gutenTAG")
#' panel_path <- system.file("extdata/ref_list.csv", package = "gutenTAG")
#' panel <- readPanel(path = panel_path)
#' raw <- readMSIData(path)
#' pre <- preProcess(raw, cores = 2)
#' peaks <- peakDetection(pre, core = 2)
#' generateMetapeaks(peaks)

generateMetapeaks <- function(x, threshold = 0.01, hist_smooth_factor = 1.0, fixed.limits = NULL) {
  count_df <- countPeaks(x)
  detection_threshold <- unname(dim(x)["Pixels"]) * threshold

  seed_mz <- generateSeedMz(smoothPeakCounts(count_df, hist_smooth_factor), detection_threshold = detection_threshold)
  metapeaks <- estimateMetapeaks(count_df, seed_mz, detection_threshold = detection_threshold, fixed.limits = fixed.limits)

  return(metapeaks)
}
