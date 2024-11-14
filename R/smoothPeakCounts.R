#' This function smoothes the peak counts histogram by convolving it with a Gaussian filter.
#' The Gaussian filter is scaled to the average step size of the mz values.
#' The function assumes that the mz values are (roughly) evenly spaced.
#'
#' @param count_df A data frame with columns: mz, count.
#' @param smooth_factor The standard deviation of the Gaussian filter is scaled by this factor.
#' @return A data frame with the same mz values and the smoothed counts.
#' @export
smoothPeakCounts <- function(count_df, smooth_factor = 1) {
  n_features <- nrow(count_df)
  sigma_smoothing_isotopic <- smooth_factor / .mzScalingFactor(count_df)
  list_values <- seq(-n_features / 2, n_features / 2, length.out = n_features)
  gaussian_vector <- exp(-list_values^2 / sigma_smoothing_isotopic^2)
  freq_peak_ordered_smooth <- stats::convolve(y = gaussian_vector, x = count_df$count, conj = TRUE)
  freq_peak_ordered_smooth <- freq_peak_ordered_smooth[c(floor(n_features / 2):n_features, 1:(n_features / 2 - 1))]

  return(data.frame(mz = count_df$mz, count = freq_peak_ordered_smooth))

  ## TODO new impl:
  # n_features <- nrow(count_df)
  # mz_vector <- count_df$mz

  ### Scale the standard deviation of the Gaussian so that it is 1 divided by the average step size.
  # sigma_smoothing_isotopic <- 1 / .mzScalingFactor(count_df)

  ### Create Gaussian filter for convolution.
  # list_values <- seq(-n_features / 2, n_features / 2, length.out = n_features)
  # gaussian_vector <- exp(-list_values^2 / sigma_smoothing_isotopic^2)

  ### Convolve the peak frequencies and the Gaussian filter
  ## TODO don't use circular, but rather "open" convolution here (breaking change)
  # conv_result <- stats::convolve(y = gaussian_vector, x = count_df$count, conj = TRUE, type = "circular")

  # .shiftN <- function(arr, n) {
  #  n_total <- length(arr)
  #  indices <- ((seq_len(n_total) + n) - 1) %% n_total + 1
  #  return(arr[indices])
  # }

  ### Shift the vector to account for the gaussian_vector being centered at index floor(n_features / 2).
  # count_smooth <- head(.shiftN(conv_result, n_features %/% 2), n_features)

  # return(data.frame(mz = mz_vector, count = count_smooth))
}
