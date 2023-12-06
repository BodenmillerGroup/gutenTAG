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
#'
#' @return A list containing information on the location, delimitations and width of metapeaks
#'
#' @import Cardinal
#' @export generateMetapeaks
#'
#' @examples
#' generateMetapeaks(list_peaks)

generateMetapeaks <- function(x, threshold = 0.01) {
  count_df <- countPeaks(x)
  detection_threshold <- unname(dim(x)["Pixels"]) * threshold

  seed_mz <- generateSeedMz(smoothPeakCounts(count_df), detection_threshold=detection_threshold)
  metapeaks <- estimateMetapeaks(count_df, seed_mz, detection_threshold=detection_threshold)
  
  return(metapeaks)
}
