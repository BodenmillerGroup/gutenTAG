#' Peak Detection
#'
#' Performs peak detection using a Median Adaptive Deviation (MAD) algorithm to discern peaks from noise.
#'
#' @param x A pre-processing MSImagingExperiment object
#' @param snr Signal to noise ratio for peak picking.
#' @param win Window size over which moving noise levels are calculated.
#' @param BPPARAM A \code{\link[BiocParallel]{BiocParallelParam}} object
#'   controlling parallel execution. Defaults to
#'   \code{\link[BiocParallel]{bpparam}()}, which uses the session-wide
#'   registered backend. Use \code{BiocParallel::MulticoreParam(n)} on
#'   Linux/macOS or \code{BiocParallel::SnowParam(n)} on Windows for
#'   multi-core processing.
#'
#' @return A list with elements \code{$peaks} (the Cardinal MSImagingExperiment
#'   object) and \code{$params} (a named list recording \code{snr} and
#'   \code{win} for downstream propagation).
#'
#' @export
#'
#' @examples
#' path <- system.file("extdata/Example_data.imzML", package = "gutenTAG")
#' panel_path <- system.file("extdata/ref_list.csv", package = "gutenTAG")
#' panel <- readPanel(path = panel_path)
#' raw <- readMSIData(path)
#' pre <- preProcess(raw)
#' peakDetection(pre)

peakDetection <- function(x, snr = 3, win = 50, BPPARAM = BiocParallel::bpparam()){

  # validity checks for peakDetection
  .valid.peakDetection(x, snr, win, BPPARAM)

  # peak detection
  picker <- Cardinal::peakPick(x, method = "mad", SNR = snr, width = win)
  list_peaks <- Cardinal::process(picker, BPPARAM = BPPARAM)

  return(list(peaks = list_peaks, params = list(snr = snr, win = win)))

}
