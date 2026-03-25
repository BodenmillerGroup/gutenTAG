#' Compute Per-Pixel Mass Shift from Preprocessed Spectral Data
#'
#' @description For each targeted metapeak, computes an intensity-weighted m/z
#'   centroid within the metapeak window for every pixel in the image. The
#'   per-pixel mass shift for marker \eqn{k} in pixel \eqn{p} is defined as
#'   the centroid m/z minus the panel reference m/z. A per-pixel consensus
#'   value is obtained by averaging the shift across all observed markers.
#'   Pixels with zero total intensity within a metapeak window receive
#'   \code{NA}.
#'
#' @param x A list. Output from \code{\link{assignMetapeaks}}, containing at
#'   minimum \code{CorrespondenceMatrix} and \code{SpatialCoords}.
#' @param pre An \code{MSImagingExperiment} object. The preprocessed spectral
#'   data returned by \code{\link{preProcess}}. Must contain the same pixels
#'   and m/z axis as the data used to generate \code{x}.
#' @param metapeaks A list. Output from \code{\link{generateMetapeaks}},
#'   providing the metapeak window limits (\code{metapeaks$metapeaks$limits})
#'   used to define the integration window per marker.
#' @param update_correspondence Logical; if \code{TRUE}, the function appends
#'   \code{MassShiftDF}, \code{PixelMassShift}, and a \code{mean_mass_shift}
#'   column to \code{x$CorrespondenceMatrix}, and returns the augmented list.
#'   If \code{FALSE} (default), returns a list with \code{MassShiftDF} and
#'   \code{PixelMassShift} directly.
#'
#' @return If \code{update_correspondence = FALSE}, a named list with two
#'   elements:
#'   \describe{
#'     \item{MassShiftDF}{Data frame of per-pixel mass shifts (pixels \eqn{\times}
#'       markers, values in Da). Columns correspond to observed (annotated)
#'       markers only.}
#'     \item{PixelMassShift}{Numeric vector of length equal to the number of
#'       pixels. Each value is the mean mass shift across all observed markers
#'       for that pixel (\code{NA} pixels are excluded via \code{na.rm = TRUE}).}
#'   }
#'   If \code{update_correspondence = TRUE}, the input list \code{x} augmented
#'   with \code{$MassShiftDF}, \code{$PixelMassShift}, and a
#'   \code{mean_mass_shift} column in \code{$CorrespondenceMatrix}.
#'
#' @examples
#' \donttest{
#'   path <- system.file("extdata/Example_data.imzML", package = "gutenTAG")
#'   panel_path <- system.file("extdata/ref_list.csv", package = "gutenTAG")
#'   panel <- readPanel(path = panel_path)
#'   raw <- readMSIData(path)
#'   pre <- preProcess(raw)
#'   peaks <- peakDetection(pre)
#'   mpeaks <- generateMetapeaks(peaks, hist_smooth_factor = 1)
#'   processed <- assignMetapeaks(mpeaks, pre = pre, refList = panel)
#'   result <- computePixelMassShift(processed, pre = pre, metapeaks = mpeaks)
#' }
#' @importFrom Cardinal mz
#' @importFrom ProtGenerics spectra
#' @importFrom methods is
#' @export computePixelMassShift
#'
computePixelMassShift <- function(x, pre, metapeaks,
                                  update_correspondence = FALSE) {

  .valid.computePixelMassShift(x, pre, metapeaks, update_correspondence)

  # Extract raw spectral data: mz_bins x pixels
  raw_intensity <- ProtGenerics::spectra(pre)
  mz_vector     <- as.numeric(Cardinal::mz(pre))
  n_pixels      <- ncol(raw_intensity)

  # Metapeak windows
  mpeaks_data <- metapeaks$metapeaks

  # Work only with observed markers (mz_location not NA)
  cm       <- x$CorrespondenceMatrix
  obs_rows <- which(!is.na(cm$mz_location))
  cm_obs   <- cm[obs_rows, ]
  n_markers <- nrow(cm_obs)

  if (n_markers == 0L) {
    stop("No observed markers found in 'x$CorrespondenceMatrix$mz_location'. ",
         "All markers are undetected; nothing to compute.")
  }

  mass_shift_mat <- matrix(NA_real_, nrow = n_pixels, ncol = n_markers)
  colnames(mass_shift_mat) <- cm_obs$marker

  for (k in seq_len(n_markers)) {

    obs_mz <- cm_obs$mz_location[k]

    # Match to metapeak by max position (within 0.01 Da tolerance)
    mp_idx <- which(abs(mpeaks_data$max - obs_mz) < 0.01)
    if (length(mp_idx) == 0L) next
    mp_idx <- mp_idx[1L]

    mz_min  <- mpeaks_data$limits[mp_idx, 1L]
    mz_max  <- mpeaks_data$limits[mp_idx, 2L]
    window  <- mz_vector >= mz_min & mz_vector <= mz_max

    if (!any(window)) next

    mz_win        <- mz_vector[window]
    intensity_win <- raw_intensity[window, , drop = FALSE]

    col_sums <- colSums(intensity_win)
    nonzero  <- col_sums > 0

    if (!any(nonzero)) next

    # Intensity-weighted centroid: sum(mz * I) / sum(I) per pixel
    weighted_sum <- as.numeric(
      matrix(mz_win, nrow = 1L) %*% intensity_win[, nonzero, drop = FALSE]
    )

    centroid_vec <- rep(NA_real_, n_pixels)
    centroid_vec[nonzero] <- weighted_sum / col_sums[nonzero]

    mass_shift_mat[, k] <- centroid_vec - cm_obs$expected_mz_location[k]
  }

  mass_shift_df    <- data.frame(mass_shift_mat, check.names = FALSE)
  pixel_mass_shift <- rowMeans(mass_shift_mat, na.rm = TRUE)
  # rowMeans returns NaN for rows that are entirely NA; normalise to NA_real_
  pixel_mass_shift[is.nan(pixel_mass_shift)] <- NA_real_

  if (update_correspondence) {
    marker_means <- colMeans(mass_shift_mat, na.rm = TRUE)
    # colMeans returns NaN for all-NA columns; normalise to NA_real_
    marker_means[is.nan(marker_means)] <- NA_real_

    x$MassShiftDF    <- mass_shift_df
    x$PixelMassShift <- pixel_mass_shift
    x$CorrespondenceMatrix$mean_mass_shift           <- NA_real_
    x$CorrespondenceMatrix$mean_mass_shift[obs_rows] <- marker_means
    return(x)
  }

  list(MassShiftDF = mass_shift_df, PixelMassShift = pixel_mass_shift)
}
