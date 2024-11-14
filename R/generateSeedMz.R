#' @title generateSeedMz
#' @description Generate seeds for the metapeak algorithm
#' @param count_df A data frame with columns mz and count
#' @param detection_threshold Number of counts below which to ignore peaks
#' @return A vector of seeds in terms of m/z values
#' @export
generateSeedMz <- function(count_df, detection_threshold = 0) {
  # constraint for how far apart metapeaks should be
  span_local_maxima <- 3 / .mzScalingFactor(count_df)

  # find peaks
  counts <- count_df$count
  counts[counts < detection_threshold] <- 0
  local_maxima <- base::which(splus2R::peaks(x = counts, span = span_local_maxima, strict = TRUE))

  # return mz values
  return(count_df$mz[local_maxima])
}
