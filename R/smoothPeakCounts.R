#' This function smoothes the peak counts histogram by convolving it with a Gaussian filter.
#' The Gaussian filter is scaled to the average step size of the mz values.
#' The function assumes that the mz values are (roughly) evenly spaced.
#'
#' @param count_df A data frame with columns: mz, count.
#' @param smooth_factor The standard deviation of the Gaussian filter is scaled by this factor.
#' @return A data frame with the same mz values and the smoothed counts.
#' @export
smoothPeakCounts <- function(count_df, smooth_factor) {

  n_features <- nrow(count_df)
  sigma_smoothing_isotopic <- smooth_factor / .mzScalingFactor(count_df)
  list_values <- seq(-n_features / 2, n_features / 2, length.out = n_features)
  gaussian_vector <- exp(-list_values^2 / sigma_smoothing_isotopic^2)
  freq_peak_ordered_smooth <- stats::convolve(y = gaussian_vector, x = count_df$count, conj = TRUE)
  freq_peak_ordered_smooth <- freq_peak_ordered_smooth[c(floor(n_features / 2):n_features, 1:(n_features/2 - 1))]

  return(data.frame(mz = count_df$mz, count = freq_peak_ordered_smooth))

  ### Convolve the peak frequencies and the Gaussian filter
  ## TODO don't use circular, but rather "open" convolution here (breaking change)
  #conv_result <- stats::convolve(y = gaussian_vector, x = count_df$count, conj = TRUE, type = "circular")

}
