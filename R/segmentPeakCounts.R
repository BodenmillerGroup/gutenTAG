#' Segment Peak Counts
#'
#' Given a data frame containing the number of peaks at a particular m/z value,
#' and a set of seed m/z values, this function segments the data frame into
#' different peaks.
#'
#' A watershed based segmentation is used.
#' A vector of class labels is returned, where each label corresponds to a peak.
#'
#' @param count_df A data frame containing counts.
#' @param seed_mz A seed mass-to-charge ratio (mz).
#' @param detection_threshold The minimum number of counts to be considered a peak.
#' @export
segmentPeakCounts <- function(count_df, seed_mz, detection_threshold) {
  # Convert the mz to row indices in the count_df
  target_indices <- findInterval(seed_mz, count_df$mz)

  # Remove values below threshold
  counts <- count_df$count
  counts[counts < detection_threshold] <- 0

  # Convert to matrix representing a column vector.
  target_indices_matrix <- matrix(0, nrow = nrow(count_df), ncol = 1)
  target_indices_matrix[target_indices] <- seq_along(target_indices)
  count_matrix <- matrix(counts, nrow = nrow(count_df), ncol = 1)

  propagation_selection <- EBImage::propagate(
    count_matrix,
    seeds = target_indices_matrix,
    mask = count_matrix > detection_threshold,
  )

  return(as.vector(propagation_selection))
}
