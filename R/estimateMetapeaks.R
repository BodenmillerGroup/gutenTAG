#' Metapeak Estimation
#'
#' This function estimates metapeaks from a given count data frame and seed points.
#' It first applies Gaussian smoothing over the counts, removes values below a certain threshold,
#' segments the histogram using seed_point m/z values, and then calculates metapeak parameters.
#' The parameters include the center, maximum, delimitation, and width of each metapeak. The width is rescaled to the m/z scale.
#'
#' @param count_df A data frame containing count data.
#' @param seed_mz A numeric vector specifying the seed m/z values for segmentation.
#' @param detection_threshold A numeric value specifying the threshold for counts, as an integer number of counts.
#' @param fixed.limits A tolerance parameter for setting the metapeak limits to be a fixed value centered around the metapeak max. The total width of the metapeak will be twice the value of this parameter.
#'
#' @return A list containing the center, maximum, delimitation, and width of each metapeak.
#'
#' @export

estimateMetapeaks <- function(count_df, seed_mz, detection_threshold = 0, fixed.limits = NULL) {
  ## validity checks
  .valid.estimateMetapeaks(count_df, seed_mz, detection_threshold)

  # 1. Gaussian smoothing over counts
  count_smooth_df <- smoothPeakCounts(count_df)

  # 2. Segment the count_df histogram using the seed_mz values.
  propagation_selection <- segmentPeakCounts(count_smooth_df, seed_mz, detection_threshold = detection_threshold)

  # 3. Get metapeak parameters
  num_metapeaks <- max(propagation_selection)
  if (num_metapeaks == 0) {
    stop("No metapeaks were found.")
  }
  metapeak_center <- double(num_metapeaks)
  metapeak_max <- double(num_metapeaks)
  metapeak_limits <- matrix(0, nrow = num_metapeaks, ncol = 2)

  for (k in seq_len(num_metapeaks)) {
    # Compute a mask of the indices that correspond to the segmentation k.
    peak_mask <- as.numeric(propagation_selection) == k
    peak_mzs <- count_df$mz[peak_mask]
    peak_counts <- count_df$count[peak_mask]

    # Compute weighted mean using peak counts.
    metapeak_center[k] <- weighted.mean(x = peak_mzs, w = peak_counts)

    # Compute maximum peak location.
    metapeak_max[k] <- peak_mzs[which.max(peak_counts)]

    if (is.null(fixed.limits)){

      # Define peak limits as being inside the 1% and 99% of the metapeak.
      cumsum_freq <- cumsum(peak_counts) / sum(peak_counts)

      # TODO more stable implementation, that does not require NA handling afterwards.
      #cumsum_rev <- rev(cumsum(rev(peak_counts)) / sum(peak_counts))
      #begin_mz <- peak_mzs[tail(which(cumsum_rev > 0.99), n = 1)[1]]
      begin_mz <- peak_mzs[tail(which(cumsum_freq < 0.01), n = 1)[1]]
      end_mz <- peak_mzs[head(which(cumsum_freq > 0.99), n = 1)[1]]

      # TODO might be worth revisiting, in principle at least the end should always be found already.
      if (!is.finite(begin_mz)) {
        begin_mz <- end_mz - 1
      }
      if (!is.finite(end_mz)) {
        end_mz <- begin_mz + 1
      }

      # Store beginning and end of peak location.
      metapeak_limits[k, ] <- c(begin_mz, end_mz)
    }else{
      # set binning limits to be within specified fixed limit around metapeak max
      metapeak_limits[k, 1] <- metapeak_max[k] - fixed.limits
      metapeak_limits[k, 2] <- metapeak_max[k] + fixed.limits
    }

  }


  # width of each peak
  # TODO is this still needed: (rescale peak width so it is on the m/z scale)
  metapeak_width <- (metapeak_limits[, 2] - metapeak_limits[, 1]) * .mzScalingFactor(count_df)
  #metapeak_width <- (metapeak_limits[, 2] - metapeak_limits[, 1])

  metapeaks <- list(center = metapeak_center,
                    max = metapeak_max,
                    width = metapeak_width,
                    limits = metapeak_limits)

  return(list(metapeaks = metapeaks,
              count_df = count_df,
              count_smooth_df = count_smooth_df,
              seed_mz = seed_mz,
              propagation_selection = propagation_selection))

}
