#' @title generateSeedMz
#' @description Generate seeds for the metapeak algorithm
#' @param count_df A data frame with columns mz and count
#' @param detection_threshold Number of counts below which to ignore peaks
#' @param sparsity Controls the minimum separation between detected peaks.
#'   Must be a single positive number. Higher values enforce greater spacing
#'   between returned seed m/z values (fewer seeds detected); lower values
#'   allow more closely-spaced seeds.
#' @return A vector of seeds in terms of m/z values
#'
#' @examples
#' rdata_path <- system.file("extdata/Example_processed.Rdata", package = "gutenTAG")
#' load(rdata_path)
#' seeds <- generateSeedMz(results$metapeaks$count_df)
#'
#' @export
generateSeedMz <- function(count_df, detection_threshold = 0, sparsity = 3) {

  if (!is.numeric(sparsity) || length(sparsity) != 1 || sparsity <= 0)
    stop("'sparsity' must be a single positive number.")

  # constraint for how far apart metapeaks should be
  span_local_maxima <- sparsity / .mzScalingFactor(count_df)

  # find peaks
  counts <- count_df$count
  counts[counts < detection_threshold] <- 0
  local_maxima <- base::which(splus2R::peaks(x = counts, span = span_local_maxima, strict = TRUE))

  # return mz values
  return(count_df$mz[local_maxima])
}
