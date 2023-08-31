#' Peak Detection
#'
#' Performs peak detection using a Median Adaptive Deviation (MAD) algorithm to discern peaks from noise.
#'
#' @param x A pre-processing MSImagingExperiment object
#' @param snr Signal to noise ratio for peak picking.
#' @param win Window size over which moving noise levels are calculated.
#' @param cores Number of cores to use for peak detection.
#'
#' @return An MSImagingExperiment object containing a list of detected peaks
#'
#' @export
#'
#' @examples
#' peakDetection(peakPre)

peakDetection <- function(x, snr = 3, win = 50, cores = 1){

  raw_intensity <- iData(x)

  range_peaks <- range(mz(x))
  n_features <- length(mz(x))
  n_pixels <- dim(x)["Pixels"]
  threshold_detection <- n_pixels*0.01
  mz_vector <- as.data.frame(mz(x))
  location_pixels <- as.data.frame(pData(x))[, c("x","y")]

  if(cores > 1){
    list_peaks <- Cardinal::peakPick(x, method = "mad", SNR = snr, window = win) %>%
      process(BPPARAM = MulticoreParam(workers = cores))
  }else{
    list_peaks <- Cardinal::peakPick(x, method = "mad", SNR = snr, window = win) %>%
      process()
  }



  return(list_peaks)

}
