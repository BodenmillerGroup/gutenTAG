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
  # Vector containing every m/z value that is present in the data, regardless of frequency of occurrence.
  # The entries are sorted and unique.
  mz_vector <- as.numeric(sort(unique(mz(mse_detected))))

  # Large list of vectors, each containing the m/z values present in a single spectrum.
  peaks_per_spectrum <- peakData(mse_detected)[['mz']]

  # Count the number of occurences of each m/z value in the list of vectors.
  peak_counts <- table(factor(unlist(as.list(peaks_per_spectrum)), levels=mz_vector))

  # Extract a data frame from the table, with mz and counts columns.
  result_df <- as.data.frame(peak_counts)
  colnames(result_df) <- c('mz', 'count')
  result_df$mz <- as.numeric(levels(result_df$mz))

  return(result_df)
}
