# Processing functions

# 1. Peak Detection ####

peakDetection <- function(x, snr = 2, win = 50){

  raw_intensity <- iData(x)

  range_peaks <- range(mz(x))
  n_features <- length(mz(x))
  n_pixels <- dim(x)["Pixels"]
  threshold_detection <- n_pixels*0.01
  mz_vector <- as.data.frame(mz(x))
  location_pixels <- as.data.frame(pData(x))[, c("x","y")]

  list_peaks <- Cardinal::peakPick(x, method = "mad", SNR = snr, window = win) %>%
    process()

  return(list_peaks)

}

test1 <- peakDetection(peakPre)




# 2. Metapeak generation

metapeakGeneration <- function(x){

  # 0. Set parameters
  mz_vector <- as.data.frame(mz(x))
  n_features <- length(mz(x))
  threshold_detection <- n_pixels*0.01

  # i. Get ordered counts
  list_peaks <- peakData(x)[["mz"]]
  n_peak_per_pixel <- unlist(lapply(list_peaks, FUN = length))

  unlist_peaks <- unlist(list_peaks)
  unlist_peaks_ordered <- unlist_peaks[order(unlist_peaks)]

  freq_peak_ordered <- table(factor(unlist_peaks_ordered, levels = mz_vector$mz))
  freq_peak_ordered <- as.numeric(freq_peak_ordered)

  # ii. Gaussian smoothing over counts
  mz_index_ratio <- (max(mz_vector$mz) - min(mz_vector$mz)) / n_features
  sigma_smoothing_isotopic <- 1 / mz_index_ratio
  list_values <- seq(-n_features / 2, n_features / 2, length.out = n_features)
  gaussian_vector <- exp(-list_values^2 / sigma_smoothing_isotopic^2)

  freq_peak_ordered_smooth <- stats::convolve(y = gaussian_vector, x = freq_peak_ordered, conj = TRUE)
  freq_peak_ordered_smooth <- freq_peak_ordered_smooth[c((n_features / 2):n_features, 1:(n_features/2 - 1))]

  # iii. Remove values below threshold
  freq_peak_smooth_thresholded <- freq_peak_ordered_smooth
  freq_peak_smooth_thresholded[freq_peak_smooth_thresholded < threshold_detection] <- 0


  # iv. Generate metapeaks
  span_local_maxima <- 3 / mz_index_ratio
  local_maxima_peak_freq <- base::which(find_peaks(freq_peak_smooth_thresholded, span = span_local_maxima))

  local_maxima_peak_freq_reshaped <- matrix(rep(0, length(freq_peak_ordered)), ncol = 1)
  local_maxima_peak_freq_reshaped[local_maxima_peak_freq] <- 1:length(local_maxima_peak_freq)

  freq_peak_smooth_thresholded <- matrix(freq_peak_smooth_thresholded, ncol = 1)

  # watershed for peak segmentation
  propagation_selection <- propagate(freq_peak_smooth_thresholded,
                                     seeds = local_maxima_peak_freq_reshaped,
                                     mask = freq_peak_smooth_thresholded > threshold_detection)


  # get metapeak location
  metapeak_location <- aggregate.data.frame(data.frame(Location = mz_vector$mz,
                                                       Freq = freq_peak_ordered),
                                            by = list(as.numeric(propagation_selection)),
                                            FUN = function(x) {sum(x[[1]]*x[[2]]) / sum(x[[2]])})
  # get metapeak parameters
  metapeak_center <- c()
  metapeak_max <- c()
  metapeak_delimitation <- c()
  for (k in 1:max(propagation_selection)) {
    # Looking for weighted mean
    location_centered_temp <- weighted.mean(x = mz_vector$mz[as.numeric(propagation_selection) == k],
                                            w = freq_peak_ordered[as.numeric(propagation_selection) == k] )
    # vector of centered peak locations
    metapeak_center <- c(metapeak_center, location_centered_temp)

    location_max_temp <- (mz_vector$mz[as.numeric(propagation_selection) == k])[which.max(freq_peak_ordered[as.numeric(propagation_selection) == k])]
    # vector of peak maximimum locations
    metapeak_max <- c(metapeak_max, location_max_temp)

    cumsum_freq <- (cumsum(freq_peak_ordered[as.numeric(propagation_selection) == k]) / sum(freq_peak_ordered[as.numeric(propagation_selection) == k]))

    beginning_peak <- max(which(cumsum_freq < 0.01))
    end_peak <- min(which(cumsum_freq > 0.99))
    if (is.infinite(beginning_peak)) {
      Beginning_peak = 1
    }
    # beginning and end of peak location
    metapeak_delimitation <- rbind(metapeak_delimitation,
                                   c(mz_vector$mz[as.numeric(propagation_selection) == k][beginning_peak],
                                     mz_vector$mz[as.numeric(propagation_selection) == k][end_peak]))
  }
  # width of each peak
  metapeak_width <- metapeak_delimitation[, 2] - metapeak_delimitation[, 1]
  # width of each peak on m/z scale
  metapeak_width <- metapeak_width * mz_index_ratio

  metapeaks <- list(center = metapeak_center,
                    max = metapeak_max,
                    width = metapeak_width,
                    limits = metapeak_delimitation)

  return(metapeaks)

}


test2 <- metapeakGeneration(test1)

getIntensityDF <- function(x, refList, pre, mz_threshold = 1){

  # get variables from pre-processed object
  mz_vector <- as.data.frame(mz(pre))
  colnames(mz_vector) <- "mz"
  raw_intensity <- iData(pre)

  mapping_meta <- crossKnn(mA = matrix(x$max),
                           mB= matrix(refList$FeatureMass, ncol = 1), k = 10, indexType = "L2", verbose = FALSE)

  # remove all mappings below the m/z distance association threshold
  mapping_meta[mapping_meta > mz_threshold] <- 0
  mapping_meta_cleaned <- apply(as.matrix(mapping_meta), MARGIN = 2, FUN = which_max_modified)

  associated_marker <- c()
  for (k in 1:length(x$max)) {
    selected_mz_location <- range(mz_vector$mz[as.numeric(propagation_selection) == k])
    annotation_assigned <- refList$Name[refList$FeatureMass > selected_mz_location[1] & refList$FeatureMass < selected_mz_location[2]]
    if (length(annotation_assigned) == 0){
      annotation_assigned = NA
    }
    if (length(annotation_assigned) > 1){
      annotation_assigned = paste(annotation_assigned, collapse="__OR__")
    }
    associated_marker <- c(associated_marker, annotation_assigned)
  }

  correspondence_matrix <- data.frame(mz_location = x$max,
                                      expected_mz_location = refList$FeatureMass[mapping_meta_cleaned],
                                      annotation_peaks = refList$Name[mapping_meta_cleaned])
  rownames(refList) <- refList$Name


  final_intensity <- c()
  for (k in 1:nrow(correspondence_matrix)) {
    # which mz locations are inside the peak (boolean)
    selected_mz_location <- mz_vector$mz >= x$limits[k,1] & mz_vector$mz <= x$limits[k,2]
    # those mz values
    mz_vector$mz[which(selected_mz_location == TRUE)]
    if (sum(selected_mz_location) > 1) {
      intensity_temp = colSums(raw_intensity[selected_mz_location, ])
    }
    # if there is only one mz location, the intensity list is only the intensity values of the one mz location
    if (sum(selected_mz_location) == 1) {
      intensity_temp <- (raw_intensity[selected_mz_location, ])
    }
    #intensity_temp = Correspondence_matrix$Annotation_peaks[k]
    final_intensity <- cbind(final_intensity, intensity_temp)

  }


  # which of the metapeaks in the final intensity matrix are annotated peaks
  final_intensity_targeted <- final_intensity[, !is.na(correspondence_matrix$Annotation_peaks)]
  colnames(final_intensity_targeted) <- correspondence_matrix$Annotation_peaks[!is.na(correspondence_matrix$Annotation_peaks)]
  final_intensity_targeted <- as.data.frame(final_intensity_targeted)

  return(correspondence_matrix, final_intensity_targeted)

}

getIntensityDF(x = test2, refList = peakAnnotation, pre = peakPre)
