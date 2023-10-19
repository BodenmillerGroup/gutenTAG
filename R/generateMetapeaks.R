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
#' @param threshold A threshold for peaks to be considered present, as a fraction of available pixels.
#'
#' @return A list containing information on the location, delimitations and width of metapeaks
#'
#' @importFrom stats convolve
#' @export generateMetapeaks
#'
#' @examples
#' generateMetapeaks(list_peaks)

generateMetapeaks <- function(x, threshold = 0.01) {

  # validity checks
  .valid.metapeakGeneration(x, threshold)

  # 0. Set parameters
  mz_vector <- as.data.frame(mz(x))
  n_pixels <- dim(x)["Pixels"]
  n_features <- length(mz(x))
  threshold_detection <- n_pixels * threshold

  # i. Count how many pixels each particular peak is present in.
  # Get a flat vector of all peaks of all pixels, including duplicates multiple times.
  unlist_peaks <- unlist(peakData(x)[["mz"]])
  # Sort the mz values of all peaks.
  unlist_peaks_ordered <- sort(unlist_peaks)

  # Count the number of instances of each of x's mz value
  freq_peak_ordered <- table(factor(unlist_peaks_ordered, levels = mz_vector$mz))
  freq_peak_ordered <- as.numeric(freq_peak_ordered)

  # ii. Gaussian smoothing over counts
  # Calculate average step size for each detected peak.
  mz_scaling_factor <- (max(mz_vector$mz) - min(mz_vector$mz)) / n_features
  # Scale the standard deviation of the Gaussian so that it is 1 divided by the average step size.
  sigma_smoothing_isotopic <- 1 / mz_scaling_factor 
  list_values <- seq(-n_features / 2, n_features / 2, length.out = n_features)
  # Create Gaussian filter for convolution.
  gaussian_vector <- exp(-list_values^2 / sigma_smoothing_isotopic^2) 

  # Convolve the peak frequencies and the Gaussian filter
  freq_peak_ordered_smooth <- stats::convolve(y = gaussian_vector, x = freq_peak_ordered, conj = TRUE)
  # Shift the vector to account for the gaussian_vector being centered at index floor(n_features / 2).
  freq_peak_ordered_smooth <- freq_peak_ordered_smooth[c(floor(n_features / 2):n_features, 1:(n_features/2 - 1))]

  # iii. Remove values below threshold
  freq_peak_smooth_thresholded <- freq_peak_ordered_smooth
  freq_peak_smooth_thresholded[freq_peak_smooth_thresholded < threshold_detection] <- 0

  # iv. Generate seeds, then metapeaks
  span_local_maxima <- 3 / mz_scaling_factor # constraint for how far apart metapeaks should be
  local_maxima_peak_freq <- base::which(.find_peaks(freq_peak_smooth_thresholded, span = span_local_maxima)) # ?

  local_maxima_peak_freq_reshaped <- matrix(rep(0, length(freq_peak_ordered)), ncol = 1)
  local_maxima_peak_freq_reshaped[local_maxima_peak_freq] <- 1:length(local_maxima_peak_freq)

  # Convert to matrix of shape (n_features, 1)
  freq_peak_smooth_thresholded <- matrix(freq_peak_smooth_thresholded, ncol = 1)

  # Segment peaks using a watershed-based algorithm.
  propagation_selection <- EBImage::propagate(freq_peak_smooth_thresholded,
                                     seeds = local_maxima_peak_freq_reshaped,
                                     mask = freq_peak_smooth_thresholded > threshold_detection)

  # Get metapeak parameters
  metapeak_center <- c()
  metapeak_max <- c()
  metapeak_delimitation <- c()

  for (k in 1:max(propagation_selection)) {
    # Compute a mask of the indices that correspond to the segmentation k.
    peak_mask <- as.numeric(propagation_selection) == k
    peak_mzs <- mz_vector$mz[peak_mask]
    peak_freqs <- freq_peak_ordered[peak_mask]

    # Compute weighted mean using peak counts.
    location_centered_temp <- weighted.mean(x = peak_mzs, w = peak_freqs)

    # Append to vector of centered peak locations.
    metapeak_center <- c(metapeak_center, location_centered_temp)

    # Append to vector of maximum peak locations.
    location_max_temp <- peak_mzs[which.max(peak_freqs)]
    metapeak_max <- c(metapeak_max, location_max_temp)

    # Define peak limits as being inside the 1% and 99% of the metapeak.
    cumsum_freq <- cumsum(peak_freqs) / sum(peak_freqs)
    begin_mz <- peak_mzs[max(which(cumsum_freq < 0.01))]
    end_mz <- peak_mzs[min(which(cumsum_freq > 0.99))]

    # If begin is NA, then compute it assuming a peak size of 1. The end should always be available by definition.
    if (is.na(begin_mz)) {
      begin_mz <- end_mz - 1
    }
    
    # Append beginning and end of peak location.
    metapeak_delimitation <- rbind(metapeak_delimitation, c(begin_mz, end_mz))
  }

  # width of each peak
  metapeak_width <- metapeak_delimitation[, 2] - metapeak_delimitation[, 1]
  # rescale peak width so it is on the m/z scale
  metapeak_width <- metapeak_width * mz_scaling_factor

  metapeaks <- list(center = metapeak_center,
                    max = metapeak_max,
                    width = metapeak_width,
                    limits = metapeak_delimitation)

  return(list(metapeaks = metapeaks,
              propagation_selection = propagation_selection))

}
