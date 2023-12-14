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
#' path <- system.file("extdata/Example_data.imzML", package = "maldipackage")
#' panel_path <- system.file("extdata/ref_list.csv", package = "maldipackage")
#' panel <- readPanel(path = panel_path)
#' raw <- readMSIData(path)
#' pre <- preProcess(raw, cores = 2)
#' peakDetection(pre)

peakDetection <- function(x, snr = 3, win = 50, cores = 1){

  # validity checks for peakDetection
  .valid.peakDetection(x, snr, win, cores)

  # is any of this necessary?
  #raw_intensity <- iData(x)

  #range_peaks <- range(mz(x))
  #n_features <- length(mz(x))
  #n_pixels <- dim(x)["Pixels"]
  #threshold_detection <- n_pixels*0.01
  #mz_vector <- as.data.frame(mz(x))
  #location_pixels <- as.data.frame(pData(x))[, c("x","y")]

  # peak detection
  if(cores > 1){
    list_peaks <- peakPick(x, method = "mad", SNR = snr, window = win) %>%
      process(BPPARAM = MulticoreParam(workers = cores))
  }else{
    list_peaks <- peakPick(x, method = "mad", SNR = snr, window = win) %>%
      process()
  }

  return(list_peaks)

}
