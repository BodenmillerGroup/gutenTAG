#' Returns the scaling factor for the mz axis.
#'
#' @param count_df A data frame with columns mz and count.
#' @return The scaling factor for the mz axis.
.mzScalingFactor <- function(count_df) {
  # TODO breaking change: divide by (nrow(count_df) + 1)
  return((max(count_df$mz) - min(count_df$mz)) / nrow(count_df))
}
