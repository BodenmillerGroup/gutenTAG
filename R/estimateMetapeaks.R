#' Metapeak Estimation
#'
#' This function estimates metapeaks from a given count data frame and seed points.
#' It first applies Gaussian smoothing over the counts, removes values below a certain threshold,
#' segments the histogram using seed_point m/z values, and then calculates metapeak parameters.
#' The parameters include the center, maximum, delimitation, and width of each metapeak. The width is rescaled to the m/z scale.
#'
#' @param count_df A data frame containing count data.
#' @param smooth_count_df A dataframe containing smooth count data.
#' @param seed_mz A numeric vector specifying the seed m/z values for segmentation.
#' @param detection_threshold A numeric value specifying the threshold for counts, as an integer number of counts.
#' @param fixed.limits A tolerance parameter for setting the metapeak limits to be a fixed value centered around the metapeak max. The total width of the metapeak will be twice the value of this parameter.
#'
#' @return A list containing the center, maximum, delimitation, and width of each metapeak.
#'
#' @examples
#' rdata_path <- system.file("extdata/Example_processed.Rdata", package = "gutenTAG")
#' load(rdata_path)
#' mp <- estimateMetapeaks(results$metapeaks$count_df,
#'                         results$metapeaks$count_smooth_df,
#'                         seed_mz = results$metapeaks$seed_mz)
#'
#' @importFrom stats weighted.mean
#' @importFrom utils head tail
#' @export

estimateMetapeaks <- function(count_df, smooth_count_df, seed_mz, detection_threshold = 0, fixed.limits = NULL) {
  ## validity checks
  .valid.estimateMetapeaks(count_df, seed_mz, detection_threshold)

  # 1. Gaussian smoothing over counts
  #count_smooth_df <- smoothPeakCounts(count_df)

  # 1. Segment the count_df histogram using the seed_mz values.
  propagation_selection <- segmentPeakCounts(smooth_count_df, seed_mz, detection_threshold = detection_threshold)

  # 2. Get metapeak parameters
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
      # Limits are the intersection of the smoothed metapeak with the
      # detection_threshold line. segmentPeakCounts() already excludes
      # bins at or below threshold, so the segment's min/max mz are
      # exactly those crossings.
      metapeak_limits[k, ] <- c(min(peak_mzs), max(peak_mzs))
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
              count_smooth_df = smooth_count_df,
              seed_mz = seed_mz,
              propagation_selection = propagation_selection))

}
