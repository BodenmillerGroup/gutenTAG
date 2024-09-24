#' Returns a sorted data.frame with a column "mz" and a column "count" indicating
#' the number of detected peaks at each mz value.
#'
#' If all counts are identical, this method will print a warning, since most likely the provided
#' MSProcessedImagingExperiment object has not been processed with peak detection.
#'
#' @param mse_detected A MSProcessedImagingExperiment object with peak detection applied to it.
#' @return data frame with mz and count columns
#'
#' @importFrom Cardinal mz peakData
#' @importFrom matter as.list
#' @export countPeaks
countPeaks <- function(mse_detected) {
  result <- Cardinal::summarizeFeatures(mse_detected, stat=c(count = "nnzero"))
  result_df <- matter::as.data.frame(Cardinal::featureData(result))
  return(result_df)
}
