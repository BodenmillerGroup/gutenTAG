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
#' rdata_path <- system.file("extdata/Example_data/Example_processed.Rdata", package = "gutenTAG")
#' load(rdata_path)
#' mp <- estimateMetapeaks(results$metapeaks$count_df,
#'                         results$metapeaks$count_smooth_df,
#'                         seed_mz = results$metapeaks$seed_mz)
#'
#' @export

estimateMetapeaks <- function(count_df, smooth_count_df, seed_mz, detection_threshold = 0, fixed.limits = NULL) {
  ## validity checks
  .valid.estimateMetapeaks(count_df, seed_mz, detection_threshold)

  # 1. Segment the count_df histogram using the seed_mz values.
  propagation_selection <- segmentPeakCounts(smooth_count_df, seed_mz, detection_threshold = detection_threshold)

  # 2. Get metapeak parameters
  num_metapeaks <- max(propagation_selection)
  if (num_metapeaks == 0) {
    stop("No metapeaks were found.")
  }
  # Each metapeak's parameters are a grouped reduction over count_df keyed by
  # propagation_selection, so compute them in one vectorised pass rather than
  # looping per metapeak. Group 0 holds the sub-threshold bins that
  # segmentPeakCounts() excluded, so drop it; pinning factor levels to
  # seq_len(num_metapeaks) keeps the output rows in metapeak order.
  g <- as.numeric(propagation_selection)
  keep <- g >= 1
  gf <- factor(g[keep], levels = seq_len(num_metapeaks))
  peak_mzs <- count_df$mz[keep]
  peak_counts <- count_df$count[keep]

  # Weighted mean m/z per metapeak: sum(mz * count) / sum(count).
  metapeak_center <- as.numeric(tapply(peak_mzs * peak_counts, gf, sum) /
                                  tapply(peak_counts, gf, sum))

  # m/z at the maximum count within each metapeak (first max, as which.max).
  metapeak_max <- as.numeric(tapply(seq_along(peak_counts), gf, function(i) {
    peak_mzs[i][which.max(peak_counts[i])]
  }))

  if (is.null(fixed.limits)) {
    # Limits are the segment's min/max m/z — the detection_threshold crossings,
    # since segmentPeakCounts() already excluded bins at or below threshold.
    metapeak_limits <- cbind(as.numeric(tapply(peak_mzs, gf, min)),
                             as.numeric(tapply(peak_mzs, gf, max)))
  } else {
    # Fixed-width window centred on each metapeak max.
    metapeak_limits <- cbind(metapeak_max - fixed.limits,
                             metapeak_max + fixed.limits)
  }


  # width of each peak (rescaled to the m/z scale)
  metapeak_width <- (metapeak_limits[, 2] - metapeak_limits[, 1]) * .mzScalingFactor(count_df)

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
