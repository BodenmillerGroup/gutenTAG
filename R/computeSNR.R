#' Compute SNR for each channel using Mixture Model (GMM).

#' @param x Output from assignMetapeaks
#' @param update_correspondence Should the SNR score be added to the Correspondence Matrix (bool)
#' @param package Indicate whether GMM modelling should be performed using mclust or flexmix package
#' @param q Quantile value for clipping extreme outliers (default is 1)
#'
#' @return a list containing lists with SNR and clustering information for each channel
#' @importFrom mclust Mclust
#' @importFrom mclust mclustBIC
#' @importFrom flexmix flexmix
#' @importFrom flexmix FLXMRglm
#' @export computeSNR
#'
#' @examples
#' path <- system.file("extdata/Example_data.imzML", package = "gutenTAG")
#' panel_path <- system.file("extdata/ref_list.csv", package = "gutenTAG")
#' panel <- readPanel(path = panel_path)
#' raw <- readMSIData(path)
#' pre <- preProcess(raw)
#' peaks <- peakDetection(pre)
#' metapeaks <- generateMetapeaks(peaks, hist_smooth_factor = 1)
#' processed <- assignMetapeaks(metapeaks, pre = pre, refList = panel)

computeSNR <- function(x, package = "mclust", update_correspondence = FALSE, q = 1) {

  # Extract data & reorder dataframe by mz order to match correspondence
  df <- x$IntensityDF[, x$CorrespondenceMatrix$marker]

  # quantile clipping in case GMM stumbles on extreme outlier pixel values
  df <- apply(df, MARGIN = 2, FUN = .quantileClipping, q = q)

  # Apply SNR calculation based on the chosen method
  if (package == "mclust") {
    sample_snrs <- apply(df, MARGIN = 2, FUN = .snr_calc_mclust)
  } else if (package == "flexmix") {
    sample_snrs <- apply(df, MARGIN = 2, FUN = .snr_calc_flexmix)
  } else {
    stop("Package must be either 'mclust' or 'flexmix'")
  }

  # Update correspondence matrix if requested
  if (update_correspondence) {
    updated <- .updateCorrespondence(x, sample_snrs)
    return(updated)
  } else {
    return(sample_snrs)
  }
}

# Helper function for updating correspondence matrix
.updateCorrespondence <- function(x, sample_snrs) {
  snrs <- vapply(sample_snrs, function(y) y$SNR, FUN.VALUE = numeric(1))
  x$CorrespondenceMatrix$snr <- snrs
  x$SNR <- sample_snrs
  return(x)
}

# Helper function: quantile clipping
.quantileClipping <- function(x, q){
  quant <- quantile(x, q)
  x[x > quant] <- quant
  return(x)
}

# SNR calculation with mclust
.snr_calc_mclust <- function(channel) {
  # Remove zero values; if channel is all zeros, return an SNR of 0
  nonzero_channel <- channel[channel != 0]
  if (length(nonzero_channel) == 0) {
    message("The channel contains only zeros. SNR set to 0.")
    return(list("SNR" = 0, "clustering" = NULL, "fit" = NULL))
  }

  # Fit a GMM with mclust's Mclust() using two components
  fit <- mclust::Mclust(log2(nonzero_channel + 1), G = 2, modelNames = "E", verbose = FALSE)  # "E" for equal variance

  if (is.null(fit)) {
    stop("The model failed to converge. Try different data or check your input.")
  }

  clusters <- fit$classification
  cluster_means <- tapply(nonzero_channel, clusters, mean)

  # Identify noise and signal clusters based on the means
  noise_cluster <- which.min(cluster_means)
  signal_cluster <- which.max(cluster_means)

  # Reassign points in signal cluster below noise mean
  threshold <- min(nonzero_channel[clusters == signal_cluster])
  clusters <- ifelse(nonzero_channel < threshold, noise_cluster, clusters)

  # Recalculate means after reassignment
  noise_mean <- mean(nonzero_channel[clusters == noise_cluster])
  signal_mean <- mean(nonzero_channel[clusters == signal_cluster])
  snr <- signal_mean / noise_mean

  return(list("SNR" = snr, "clustering" = clusters, "fit" = fit))
}

# SNR calculation with flexmix
.snr_calc_flexmix <- function(channel) {
  nonzero_channel <- channel[channel != 0]
  if (length(nonzero_channel) == 0) {
    message("The channel contains only zeros. SNR set to 0.")
    return(list("SNR" = 0, "clustering" = NULL, "fit" = NULL))
  }

  data <- data.frame(nonzero_channel)
  fit <- tryCatch({
    flexmix(log2(nonzero_channel + 1) ~ 1, data = data, k = 2,
            model = FLXMRglm(family = "gaussian"),
            control = list(iter.max = 500, tol = 1e-6))
  }, error = function(e) {
    stop("flexmix model failed to converge: ", conditionMessage(e), call. = FALSE)
  })

  clusters <- flexmix::clusters(fit)
  cluster_means <- tapply(nonzero_channel, clusters, mean)

  # Identify noise and signal clusters based on the means
  noise_cluster <- which.min(cluster_means)
  signal_cluster <- which.max(cluster_means)

  # Reassign points in signal cluster below noise mean
  threshold <- min(nonzero_channel[clusters == signal_cluster])
  clusters <- ifelse(nonzero_channel < threshold, noise_cluster, clusters)

  # Recalculate means after reassignment
  noise_mean <- mean(nonzero_channel[clusters == noise_cluster])
  signal_mean <- mean(nonzero_channel[clusters == signal_cluster])
  snr <- signal_mean / noise_mean

  return(list("SNR" = snr, "clustering" = clusters, "fit" = fit))
}

