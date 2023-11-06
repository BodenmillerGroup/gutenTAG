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
#'
#' @return A list containing information on the location, delimitations and width of metapeaks
#'
#' @importFrom stats convolve
#' @importFrom ggpmisc find_peaks
#' @export generateMetapeaks
#'
#' @examples
#' generateMetapeaks(list_peaks)

generateMetapeaks <- function(x, threshold = 0.01){

  # validity checks
  .valid.metapeakGeneration(x, threshold)

  # 0. Set parameters
  mz_vector <- as.data.frame(mz(x))
  n_pixels <- dim(x)["Pixels"]
  n_features <- length(mz(x))
  threshold_detection <- n_pixels * threshold

  # i. Get ordered counts
  list_peaks <- peakData(x)[["mz"]]
  n_peak_per_pixel <- unlist(lapply(list_peaks, FUN = length))

  unlist_peaks <- unlist(list_peaks)
  unlist_peaks_ordered <- unlist_peaks[order(unlist_peaks)]

  freq_peak_ordered <- table(factor(unlist_peaks_ordered, levels = mz_vector$mz))
  freq_peak_ordered <- as.numeric(freq_peak_ordered)

  # ii. Gaussian smoothing over counts
  mz_scaling_factor <- (max(mz_vector$mz) - min(mz_vector$mz)) / n_features # calculate avg step size for each detected peak
  sigma_smoothing_isotopic <- 1 / mz_scaling_factor # scale the standard deviation of the Gaussian so that it is 1 divided by the avg step size
  list_values <- seq(-n_features / 2, n_features / 2, length.out = n_features)
  gaussian_vector <- exp(-list_values^2 / sigma_smoothing_isotopic^2) # create Gaussian filter for convolution

  # convolve the peak frequencies and the Gaussian filter
  freq_peak_ordered_smooth <- stats::convolve(y = gaussian_vector, x = freq_peak_ordered, conj = TRUE)
  freq_peak_ordered_smooth <- freq_peak_ordered_smooth[c((n_features / 2):n_features, 1:(n_features/2 - 1))]

  # iii. Remove values below threshold
  freq_peak_smooth_thresholded <- freq_peak_ordered_smooth
  freq_peak_smooth_thresholded[freq_peak_smooth_thresholded < threshold_detection] <- 0


  # iv. Generate seeds, then metapeaks
  span_local_maxima <- 3 / mz_scaling_factor # constraint for how far apart metapeaks should be
  local_maxima_peak_freq <- base::which(.find_peaks(freq_peak_smooth_thresholded, span = span_local_maxima)) #

  #local_maxima_peak_freq <- base::which(ggpmisc:::find_peaks(freq_peak_smooth_thresholded, span = span_local_maxima)) # find_peaks from ggpmisc package



  local_maxima_peak_freq_reshaped <- matrix(rep(0, length(freq_peak_ordered)), ncol = 1)
  local_maxima_peak_freq_reshaped[local_maxima_peak_freq] <- 1:length(local_maxima_peak_freq)

  freq_peak_smooth_thresholded <- matrix(freq_peak_smooth_thresholded, ncol = 1)

  # watershed for peak segmentation
  propagation_selection <- EBImage::propagate(freq_peak_smooth_thresholded,
                                     seeds = local_maxima_peak_freq_reshaped,
                                     mask = freq_peak_smooth_thresholded > threshold_detection)


  # get metapeak location
  metapeak_location <- aggregate.data.frame(data.frame(Location = mz_vector$mz,
                                                       Freq = freq_peak_ordered),
                                            by = list(as.numeric(propagation_selection)),
                                            FUN = function(x) {sum(x[[1]]*x[[2]]) / sum(x[[2]])})
  # get metapeak parameters
  metapeak_center <- c()
  metapeak_max <- c()
  metapeak_delimitation <- c()
  for (k in 1:max(propagation_selection)) {
    # Looking for weighted mean
    location_centered_temp <- weighted.mean(x = mz_vector$mz[as.numeric(propagation_selection) == k],
                                            w = freq_peak_ordered[as.numeric(propagation_selection) == k] )
    # vector of centered peak locations
    metapeak_center <- c(metapeak_center, location_centered_temp)

    location_max_temp <- (mz_vector$mz[as.numeric(propagation_selection) == k])[which.max(freq_peak_ordered[as.numeric(propagation_selection) == k])]
    # vector of peak maximimum locations
    metapeak_max <- c(metapeak_max, location_max_temp)

    cumsum_freq <- (cumsum(freq_peak_ordered[as.numeric(propagation_selection) == k]) / sum(freq_peak_ordered[as.numeric(propagation_selection) == k]))

    # define peak limits as being inside the 1% and 99% of the metapeak
    beginning_peak <- max(which(cumsum_freq < 0.01))
    end_peak <- min(which(cumsum_freq > 0.99)) #
    if (is.infinite(beginning_peak)) {
      Beginning_peak = 1
    }
    # beginning and end of peak location
    metapeak_delimitation <- rbind(metapeak_delimitation,
                                   c(mz_vector$mz[as.numeric(propagation_selection) == k][beginning_peak],
                                     mz_vector$mz[as.numeric(propagation_selection) == k][end_peak]))
  }

  # replace NA values in metapeak_delimitation with values 1mz lower or higher than the existing value for the other limit
  ## creates an artificial peak width of 1mz (before m/z scaling) if it can not be directly calculated
  for (i in 1:nrow(metapeak_delimitation)){

    if(is.na(metapeak_delimitation[i, 1])){
      metapeak_delimitation[i, 1] <- metapeak_delimitation[i, 2] - 1 # replace peak-beginning NA with peak-end minus 1
    }

    if(is.na(metapeak_delimitation[i, 2])){
      metapeak_delimitation[i, 2] <- metapeak_delimitation[i, 1] + 1 # replace peak-end NA with peak-beginning NA plus 1
    }

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
