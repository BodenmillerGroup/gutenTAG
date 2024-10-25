#' Compute SNR for each channel using Mixture Model (GMM).

#' @param x Output from assignMetapeaks
#' @param method Indicate whether to use Gaussian mixture model or Poisson mixture model
#' @param update_correspondence Should the SNR score be added to the Correspondence Matrix (bool)
#'
#' @return a list containing lists with SNR and clustering information for each channel
#' @export computeSNR
#'
#' @examples
#' path <- system.file("extdata/Example_data.imzML", package = "gutenTAG")
#' panel_path <- system.file("extdata/ref_list.csv", package = "gutenTAG")
#' panel <- readPanel(path = panel_path)
#' raw <- readMSIData(path)
#' pre <- preProcess(raw, cores = 2)
#' peaks <- peakDetection(pre, core = 2)
#' metapeaks <- generateMetapeaks(peaks)
#' processed <- assignMetapeaks(metapeaks, pre = pre, refList = panel)

computeSNR <- function(x, method = "gaussian", update_correspondence = FALSE){

  # extract data & reorder dataframe by mz order to match correspondence
  df <- x$IntensityDF[, x$CorrespondenceMatrix$marker]

  # apply SNR to dataframe
  sample_snrs <- apply(df, MARGIN = 2, FUN = .snr_calc)

  .updateCorrespondence <- function(x, sample_snrs){
    # extract snr values only
    snrs <- sapply(sample_snrs, function(y) y$SNR)
    # add snrs to correspondence matrix snr column
    x$CorrespondenceMatrix$snr <- snrs
    # add all snr data to processed object
    x$SNR <- sample_snrs

    return(x)
  }

  # update correspondence or not
  if(update_correspondence == TRUE){
    updated <- .updateCorrespondence(x, sample_snrs)
    return(updated)
  }else{
    return(sample_snrs)
  }

}

# hidden functions
.snr_calc <- function(channel, method = "gaussian") {

  # Remove zero values, if channel is all zeros insert dummy
  nonzero_channel <- channel[channel != 0]

  if (length(nonzero_channel) == 0) {
    message("The channel contains only zeros. SNR set to 0.")
    return(list("SNR" = 0, "clustering" = NULL, "fit" = NULL))
  }

  # Convert to a data frame
  data <- data.frame(nonzero_channel)

  # Try Gaussian or Poisson mixture
  if (method == "gaussian") {
    fit <- tryCatch({
      flexmix::flexmix(log2(nonzero_channel + 1) ~ 1, data = data, k = 2,
              model = flexmix::FLXMRglm(family = "gaussian"),
              control = list(iter.max = 500, tol = 1e-6))
    }, error = function(e) {
      message("Error in Gaussian fitting: ", e)
      return(NULL)
    })

  } else if (method == "poisson") {
    fit <- tryCatch({
      flexmix::flexmix(nonzero_channel ~ 1, data = data, k = 2,
              model = flexmix::FLXMRglm(family = "poisson"),
              control = list(iter.max = 500, tol = 1e-6))
    }, error = function(e) {
      message("Error in Poisson fitting: ", e)
      return(NULL)
    })

  } else {
    stop("Method must be either 'gaussian' or 'poisson'.")
  }

  # check if fitting worked
  if (is.null(fit)) {
    stop("The model failed to converge. Try different methods or check your data.")
  }

  # Get signal and noise clusters
  clusters <- flexmix::clusters(fit)

  # Calculate the mean for each signal and noise clusters
  cluster_means <- tapply(nonzero_channel, clusters, mean)
  noise_mean <- min(cluster_means)
  signal_mean <- max(cluster_means)

  # Compute SNR (Signal-to-Noise Ratio)
  snr <- signal_mean / noise_mean

  out <- list("SNR" = snr,
              "clustering" = clusters,
              "fit" = fit
              )

  return(out)
}

