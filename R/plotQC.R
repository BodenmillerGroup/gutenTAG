# QC Plot Functions for gutenTAG
# All functions in this file produce quality control plots for MALDI-MSI data.

#' @importFrom ggplot2 ggplot aes geom_bar geom_boxplot geom_histogram
#'   geom_hline geom_line geom_point geom_rect geom_violin geom_vline
#'   scale_color_manual scale_fill_identity scale_fill_manual scale_linetype_manual
#'   scale_x_log10 scale_y_log10 labs ggtitle theme theme_minimal theme_void
#'   element_text
#' @importFrom rlang .data
#' @importFrom dplyr arrange filter
#' @importFrom stats lm sd quantile
NULL


# ── 1. plotMetapeaks ──────────────────────────────────────────────────────────

#' Plot Metapeak Segmentation
#'
#' @description Dual-axis plot showing the peak count histogram (line), smoothed
#'   counts, metapeak limit rectangles, panel feature masses (vertical lines),
#'   and optionally mean and skyline spectra overlaid on a secondary axis
#'   (interactive mode only, via \pkg{plotly}).
#'
#' @param x List output of \code{assignMetapeaks}.
#' @param metapeaks List output of \code{generateMetapeaks}.
#' @param panel Data frame from \code{readPanel} with columns \code{FeatureMass}
#'   and \code{Name}.
#' @param interactive Logical; if \code{TRUE} returns a \code{plotly} object
#'   with mean and skyline spectra on a secondary y-axis. Requires the
#'   \pkg{plotly} package. Default \code{FALSE}.
#'
#' @return A \code{ggplot} object, or a \code{plotly} object when
#'   \code{interactive = TRUE}.
#'
#' @examples
#' \dontrun{
#'   plotMetapeaks(x = processed, metapeaks = metapeaks, panel = panel)
#'   plotMetapeaks(x = processed, metapeaks = metapeaks, panel = panel,
#'                 interactive = TRUE)
#' }
#' @export
plotMetapeaks <- function(x, metapeaks, panel, interactive = FALSE) {

  # Build metapeak data frame
  mp <- metapeaks$metapeaks
  metapeak_df <- data.frame(
    center      = mp$center,
    max         = mp$max,
    lower_limit = mp$limits[, 1],
    upper_limit = mp$limits[, 2]
  )

  max_y <- max(metapeaks$count_df$count, na.rm = TRUE)

  # Identify targeted metapeaks
  targeted_mz <- x$CorrespondenceMatrix$mz_location
  targeted_metapeak_df <- metapeak_df[metapeak_df$max %in% targeted_mz, ]

  p <- ggplot() +
    theme_minimal() +
    ggtitle("Metapeaks") +
    geom_line(
      data = metapeaks$count_df,
      aes(x = .data[["mz"]], y = .data[["count"]], color = "Count"),
      alpha = 1, linewidth = 0.5
    ) +
    geom_line(
      data = metapeaks$count_smooth_df,
      aes(x = .data[["mz"]], y = .data[["count"]] / 8, color = "SmoothCount"),
      alpha = 1, linewidth = 0.5
    ) +
    geom_vline(
      data = panel,
      aes(xintercept = .data[["FeatureMass"]], linetype = "Feature Mass"),
      color = "purple", linewidth = 1
    ) +
    geom_vline(
      data = metapeak_df,
      aes(xintercept = .data[["max"]], linetype = "Max Counts"),
      color = "springgreen3", linewidth = 1
    ) +
    geom_vline(
      aes(xintercept = metapeaks$seed_mz, linetype = "Seeds"),
      color = "turquoise", linewidth = 1
    ) +
    geom_rect(
      data = metapeak_df,
      aes(xmin = .data[["lower_limit"]], xmax = .data[["upper_limit"]],
          ymin = 0, ymax = max_y, fill = "Metapeak Range"),
      alpha = 0.2
    ) +
    geom_rect(
      data = targeted_metapeak_df,
      aes(xmin = .data[["lower_limit"]], xmax = .data[["upper_limit"]],
          ymin = 0, ymax = max_y, fill = "Targeted Metapeak Range"),
      alpha = 0.2
    ) +
    labs(y = "Peak Counts") +
    scale_color_manual(
      name   = "Line",
      values = c("Count" = "black", "SmoothCount" = "sienna")
    ) +
    scale_linetype_manual(
      name   = "Vertical Lines",
      values = c("Feature Mass" = "solid", "Max Counts" = "dashed", "Seeds" = "dashed")
    ) +
    scale_fill_manual(
      name   = "Range Areas",
      values = c("Metapeak Range" = "grey", "Targeted Metapeak Range" = "red3")
    )

  if (!interactive) {
    return(p)
  }

  if (!requireNamespace("plotly", quietly = TRUE)) {
    stop(
      "Package 'plotly' is required for interactive = TRUE. ",
      "Install it with install.packages('plotly')."
    )
  }

  plotly_p <- plotly::ggplotly(p, dynamicTicks = TRUE)

  plotly_p <- plotly::add_trace(
    plotly_p,
    data = x$SummarySpectra,
    x = ~mz, y = ~mean,
    type = "scatter", mode = "lines",
    name = "MeanSpectrum",
    yaxis = "y2",
    line = list(color = "#cd0000", width = 2)
  )
  plotly_p <- plotly::add_trace(
    plotly_p,
    data = x$SummarySpectra,
    x = ~mz, y = ~skyline,
    type = "scatter", mode = "lines",
    name = "SkylineSpectrum",
    yaxis = "y2",
    line = list(color = "blue", width = 2)
  )
  plotly_p <- plotly::layout(
    plotly_p,
    yaxis  = list(
      title     = "Peak Counts",
      titlefont = list(size = 17, color = "black"),
      tickfont  = list(size = 14, color = "black")
    ),
    yaxis2 = list(
      title      = "Spectral Intensity",
      overlaying = "y",
      side       = "right",
      showgrid   = FALSE,
      titlefont  = list(size = 17, color = "#cd0000"),
      tickfont   = list(size = 14, color = "#cd0000"),
      tickpadding = 0
    ),
    margin = list(r = 1),
    legend = list(
      x = 1.15, y = 1,
      title      = list(text = "Legend"),
      showlegend = TRUE,
      font       = list(size = 8)
    )
  )

  return(plotly_p)
}


# ── 2. plotIntensityDistribution ──────────────────────────────────────────────

#' Plot Intensity Distribution per Marker
#'
#' @description Violin and boxplot showing the distribution of intensity values
#'   per marker, ordered by descending mean intensity.
#'
#' @param x List output of \code{assignMetapeaks}.
#'
#' @return A \code{ggplot} object.
#'
#' @examples
#' \dontrun{
#'   plotIntensityDistribution(x = processed)
#' }
#' @export
plotIntensityDistribution <- function(x) {

  df <- x$IntensityDF
  df_long <- utils::stack(df)
  colnames(df_long) <- c("intensity", "marker")

  ggplot(df_long, aes(x = reorder(.data[["marker"]], -.data[["intensity"]]),
                      y = .data[["intensity"]])) +
    geom_violin(
      fill  = "grey80",
      color = "black",
      scale = "width",
      trim  = TRUE
    ) +
    geom_boxplot(
      width         = 0.15,
      outlier.shape = NA,
      fill          = "white",
      color         = "black"
    ) +
    theme_minimal() +
    labs(
      x     = "",
      y     = "Intensity per pixel",
      title = "Mean Marker Intensities"
    ) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
}


# ── 3. plotMeanVariance ───────────────────────────────────────────────────────

#' Plot Mean-Variance Relationship
#'
#' @description Log-log scatter of mean vs standard deviation for all
#'   metapeaks. Targeted metapeaks are highlighted in red. When
#'   \code{annotated_only = TRUE}, only targeted metapeaks are shown in blue.
#'
#' @param x List output of \code{assignMetapeaks}.
#' @param annotated_only Logical; if \code{TRUE}, plot only targeted metapeaks.
#'   Default \code{FALSE}.
#' @param interactive Logical; if \code{TRUE} returns a \code{plotly} object.
#'   Requires the \pkg{plotly} package. Default \code{FALSE}.
#'
#' @return A \code{ggplot} object, or a \code{plotly} object when
#'   \code{interactive = TRUE}.
#'
#' @examples
#' \dontrun{
#'   plotMeanVariance(x = processed)
#'   plotMeanVariance(x = processed, annotated_only = TRUE, interactive = TRUE)
#' }
#' @export
plotMeanVariance <- function(x, annotated_only = FALSE, interactive = FALSE) {

  untargeted_df   <- data.frame(x$Untargeted$UntargetedIntensity)
  untargeted_mean <- colMeans(untargeted_df)
  untargeted_sd   <- apply(untargeted_df, MARGIN = 2, FUN = sd)
  untargeted_stats <- data.frame(
    untargeted_mean = untargeted_mean,
    untargeted_sd   = untargeted_sd
  )

  targeted_mean <- colMeans(x$IntensityDF)
  targeted_sd   <- apply(x$IntensityDF, MARGIN = 2, FUN = sd)
  targeted_stats <- data.frame(
    targeted_mean = targeted_mean,
    targeted_sd   = targeted_sd,
    Name          = rownames(data.frame(targeted_mean))
  )

  if (annotated_only) {
    p <- ggplot() +
      theme_minimal() +
      ggtitle("Mean-Variance Plot: Annotated Peaks Only") +
      geom_point(
        data  = targeted_stats,
        aes(x = .data[["targeted_mean"]], y = .data[["targeted_sd"]],
            text = .data[["Name"]]),
        color = "black", fill = "#3B9AB5", shape = 21, size = 4
      ) +
      scale_x_log10() +
      scale_y_log10() +
      labs(x = "Mean", y = "SD")
  } else {
    p <- ggplot() +
      theme_minimal() +
      ggtitle("Mean-Variance Plot") +
      geom_point(
        data  = untargeted_stats,
        aes(x = .data[["untargeted_mean"]], y = .data[["untargeted_sd"]]),
        color = "black", fill = "grey", shape = 21, size = 2
      ) +
      geom_point(
        data  = targeted_stats,
        aes(x = .data[["targeted_mean"]], y = .data[["targeted_sd"]],
            text = paste0("Marker: ", .data[["Name"]])),
        color = "black", fill = "red3", shape = 21, size = 3
      ) +
      scale_x_log10() +
      scale_y_log10() +
      labs(x = "Mean", y = "SD")
  }

  if (!interactive) {
    return(p)
  }

  if (!requireNamespace("plotly", quietly = TRUE)) {
    stop(
      "Package 'plotly' is required for interactive = TRUE. ",
      "Install it with install.packages('plotly')."
    )
  }

  plotly::ggplotly(p, dynamicTicks = TRUE)
}


# ── 4. plotMeanVarianceResiduals ───────────────────────────────────────────────

#' Plot Mean-Variance Residuals
#'
#' @description Scatter plot of residuals from a linear fit of
#'   \code{log(1 + sd)} against \code{log(1 + mean)} across all metapeaks.
#'   A horizontal reference line at zero is drawn in red. Annotated
#'   (targeted) metapeaks are coloured \code{red3}; untargeted metapeaks
#'   are shown in grey. Useful for identifying outlier markers that deviate
#'   from the mean-variance trend.
#'
#' @param x List output of \code{\link{assignMetapeaks}}.
#' @param standardised Logical; if \code{TRUE}, the y-axis shows
#'   \code{residuals / sd} (dimensionless standardised residuals).
#'   Default \code{FALSE}.
#' @param interactive Logical; if \code{TRUE} returns a \code{plotly} object.
#'   Requires the \pkg{plotly} package. Default \code{FALSE}.
#'
#' @return A \code{ggplot} object, or a \code{plotly} object when
#'   \code{interactive = TRUE}.
#'
#' @examples
#' \dontrun{
#'   plotMeanVarianceResiduals(x = processed)
#'   plotMeanVarianceResiduals(x = processed, standardised = TRUE)
#'   plotMeanVarianceResiduals(x = processed, interactive = TRUE)
#' }
#' @export
plotMeanVarianceResiduals <- function(x, standardised = FALSE,
                                      interactive = FALSE) {

  all_df         <- data.frame(x$Untargeted$UntargetedIntensity,
                                check.names = FALSE)
  mean_intensity <- colMeans(all_df)
  sd_intensity   <- apply(all_df, MARGIN = 2, FUN = sd)

  lm_fit    <- lm(log(1 + sd_intensity) ~ log(1 + mean_intensity))
  raw_resid <- lm_fit$residuals
  std_resid <- raw_resid / sd_intensity

  resid_df <- data.frame(
    mean      = mean_intensity,
    residual  = if (standardised) std_resid else raw_resid,
    annotated = colnames(all_df) %in%
                  paste0(round(x$CorrespondenceMatrix$mz_location, 2), " m/z"),
    label     = colnames(all_df),
    stringsAsFactors = FALSE
  )

  annotated_idx <- which(resid_df$annotated)
  cm_labels <- setNames(x$CorrespondenceMatrix$marker,
                        paste0(round(x$CorrespondenceMatrix$mz_location, 2), " m/z"))
  resid_df$label[annotated_idx] <- cm_labels[resid_df$label[annotated_idx]]

  y_label <- if (standardised) "Standardised residuals" else "Residuals"
  title   <- if (standardised) "Standardised Residuals vs Mean Intensity" else
               "Residuals vs Mean Intensity"

  untargeted_df <- resid_df[!resid_df$annotated, ]
  targeted_df   <- resid_df[ resid_df$annotated, ]

  p <- ggplot() +
    theme_minimal() +
    ggtitle(title) +
    geom_hline(yintercept = 0, color = "red", linetype = "solid") +
    geom_point(
      data  = untargeted_df,
      aes(x = log(1 + .data[["mean"]]), y = .data[["residual"]]),
      color = "black", fill = "grey", shape = 21, size = 2
    ) +
    geom_point(
      data  = targeted_df,
      aes(x = log(1 + .data[["mean"]]),
          y = .data[["residual"]],
          text = paste0("Marker: ", .data[["label"]])),
      color = "black", fill = "red3", shape = 21, size = 3
    ) +
    labs(x = "log(1 + mean)", y = y_label)

  if (!interactive) return(p)

  if (!requireNamespace("plotly", quietly = TRUE)) {
    stop(
      "Package 'plotly' is required for interactive = TRUE. ",
      "Install it with install.packages('plotly')."
    )
  }

  plotly::ggplotly(p, dynamicTicks = TRUE)
}


# ── 5. plotGearysC ────────────────────────────────────────────────────────────

#' Plot Geary's C Score per Marker
#'
#' @description Bar chart of Geary's C spatial autocorrelation score per marker,
#'   sorted ascending, with a dashed horizontal line at the mean score.
#'
#' @param x List output of \code{computeGearysC(..., update_correspondence = TRUE)}.
#'   \code{x$CorrespondenceMatrix} must contain a \code{GearysC} column.
#'
#' @return A \code{ggplot} object.
#'
#' @examples
#' \dontrun{
#'   processed_geary <- computeGearysC(processed, update_correspondence = TRUE)
#'   plotGearysC(x = processed_geary)
#' }
#' @export
plotGearysC <- function(x) {

  cm <- x$CorrespondenceMatrix
  plot_df <- cm[!is.na(cm$GearysC), ]
  plot_df <- plot_df[order(plot_df$GearysC, decreasing = FALSE), ]
  geary_mean <- mean(plot_df$GearysC)

  ggplot(plot_df, aes(x = reorder(.data[["marker"]], .data[["GearysC"]]),
                      y = .data[["GearysC"]])) +
    geom_bar(stat = "identity", fill = "#b51837", color = "black", linewidth = 0.2) +
    geom_hline(yintercept = geary_mean, linetype = "dashed", color = "grey") +
    theme_minimal() +
    labs(
      title = "Geary's C Score",
      x     = NULL,
      y     = "Geary's C"
    ) +
    theme(
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 7)
    )
}


# ── 6. plotSNR ────────────────────────────────────────────────────────────────

#' Plot Signal-to-Noise Ratio per Marker
#'
#' @description Bar chart of SNR per marker sorted ascending. Bars are coloured
#'   by whether their SNR falls below the threshold (red) or above it (teal).
#'
#' @param x List output of \code{computeSNR(..., update_correspondence = TRUE)}.
#' @param snr_threshold Numeric; threshold for fill colour. If \code{NULL}
#'   (default), uses the 10th percentile of SNR values.
#'
#' @return A \code{ggplot} object.
#'
#' @examples
#' \dontrun{
#'   processed_snr <- computeSNR(processed, update_correspondence = TRUE)
#'   plotSNR(x = processed_snr)
#'   plotSNR(x = processed_snr, snr_threshold = 5)
#' }
#' @export
plotSNR <- function(x, snr_threshold = NULL) {

  cm <- x$CorrespondenceMatrix
  plot_df <- cm[!is.na(cm$mz_location), ]
  plot_df <- dplyr::arrange(plot_df, .data[["snr"]])

  if (is.null(snr_threshold)) {
    snr_threshold <- quantile(plot_df$snr, 0.1)
  }

  ggplot(plot_df,
         aes(x = reorder(.data[["marker"]], -.data[["snr"]]),
             y = .data[["snr"]],
             fill = .data[["snr"]] <= snr_threshold)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(
      values = c("FALSE" = "#18B596", "TRUE" = "#b51837"),
      labels = c("FALSE" = "Above threshold", "TRUE" = "Below threshold"),
      name   = NULL
    ) +
    geom_hline(yintercept = snr_threshold, linetype = "dashed", color = "black") +
    labs(
      title = "SNR per Marker",
      x     = NULL,
      y     = "SNR"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 7)
    )
}


# ── 6. plotSNRHistogram ───────────────────────────────────────────────────────

#' Plot SNR Histogram for a Single Channel
#'
#' @description Histogram of intensity values for a single channel, with bars
#'   coloured by signal (red) vs noise (grey) classification from
#'   \code{computeSNR}.
#'
#' @param x List output of \code{computeSNR(..., update_correspondence = TRUE)}.
#' @param channel Integer index or character name of the channel to plot.
#'
#' @return A \code{ggplot} object.
#'
#' @examples
#' \dontrun{
#'   processed_snr <- computeSNR(processed, update_correspondence = TRUE)
#'   plotSNRHistogram(x = processed_snr, channel = 1)
#'   plotSNRHistogram(x = processed_snr, channel = "CD4")
#' }
#' @export
plotSNRHistogram <- function(x, channel) {

  chan_values <- x$IntensityDF[, channel]

  if (is.numeric(channel)) {
    chan_name <- colnames(x$IntensityDF)[channel]
  } else {
    chan_name <- channel
  }

  nonzero_channel <- chan_values[chan_values != 0]
  quant <- quantile(nonzero_channel, 0.999)
  nonzero_channel[nonzero_channel > quant] <- quant

  cur_cluster <- x$SNR[[chan_name]]$clustering
  cur_snr     <- x$SNR[[chan_name]]$SNR

  cluster_means <- tapply(nonzero_channel, cur_cluster, mean)
  cluster_color <- ifelse(
    cur_cluster == which.min(cluster_means), "Noise", "Signal"
  )

  hist_df <- data.frame(
    count         = nonzero_channel,
    cluster_color = cluster_color
  )

  ggplot(hist_df, aes(x = .data[["count"]], fill = .data[["cluster_color"]])) +
    geom_histogram(bins = 100, color = "black", alpha = 0.5,
                   position = "identity") +
    scale_fill_manual(
      values = c("Noise" = "grey", "Signal" = "red"),
      name   = "Cluster"
    ) +
    theme_minimal() +
    labs(
      title    = paste("Histogram of Signal and Noise:", chan_name),
      subtitle = paste("Signal-to-Noise Ratio (SNR):", round(cur_snr, 2)),
      x        = "Intensity",
      y        = "Count"
    )
}


# ── 7. plotSpatialSNR ─────────────────────────────────────────────────────────

#' Plot Spatial Distribution of Signal and Noise for a Single Channel
#'
#' @description Spatial scatter plot (x/y pixel coordinates) coloured by signal
#'   vs noise classification from \code{computeSNR} for a single channel.
#'
#' @param x List output of \code{computeSNR(..., update_correspondence = TRUE)}.
#' @param channel Integer index or character name of the channel to plot.
#'
#' @return A \code{ggplot} object.
#'
#' @examples
#' \dontrun{
#'   processed_snr <- computeSNR(processed, update_correspondence = TRUE)
#'   plotSpatialSNR(x = processed_snr, channel = 1)
#'   plotSpatialSNR(x = processed_snr, channel = "CD4")
#' }
#' @export
plotSpatialSNR <- function(x, channel) {

  chan_values <- x$IntensityDF[, channel]

  if (is.numeric(channel)) {
    chan_name <- colnames(x$IntensityDF)[channel]
  } else {
    chan_name <- channel
  }

  nonzero_mask    <- chan_values != 0
  nonzero_channel <- chan_values[nonzero_mask]
  coords          <- x$SpatialCoords[nonzero_mask, ]

  cur_cluster <- x$SNR[[chan_name]]$clustering
  cur_snr     <- x$SNR[[chan_name]]$SNR

  cluster_means <- tapply(nonzero_channel, cur_cluster, mean)
  cluster_color <- ifelse(
    cur_cluster == which.min(cluster_means), "Noise", "Signal"
  )

  plot_df <- data.frame(coords, cluster_color = cluster_color)

  ggplot(plot_df,
         aes(x = .data[["x"]], y = -.data[["y"]],
             color = factor(.data[["cluster_color"]]))) +
    geom_point(size = 7) +
    scale_color_manual(
      values = c("Noise" = "grey", "Signal" = "red"),
      name   = "Cluster"
    ) +
    theme_void() +
    labs(
      title = paste("Spatial Distribution:", chan_name),
      x     = "X",
      y     = "Y"
    )
}


# ── 8. plotMeanVsSNR ──────────────────────────────────────────────────────────

#' Plot Mean Intensity vs SNR
#'
#' @description Scatter plot of per-marker mean intensity vs signal-to-noise
#'   ratio (SNR).
#'
#' @param x List output of \code{computeSNR(..., update_correspondence = TRUE)}.
#'
#' @return A \code{ggplot} object.
#'
#' @examples
#' \dontrun{
#'   processed_snr <- computeSNR(processed, update_correspondence = TRUE)
#'   plotMeanVsSNR(x = processed_snr)
#' }
#' @export
plotMeanVsSNR <- function(x) {

  targeted_mean <- colMeans(x$IntensityDF)
  targeted_snr  <- x$CorrespondenceMatrix$snr
  plot_df <- data.frame(
    targeted_mean = targeted_mean,
    targeted_snr  = targeted_snr,
    Name          = colnames(x$IntensityDF)
  )

  ggplot() +
    theme_minimal() +
    ggtitle("Mean-SNR Plot") +
    geom_point(
      data  = plot_df,
      aes(x    = .data[["targeted_mean"]],
          y    = .data[["targeted_snr"]],
          text = .data[["Name"]]),
      color = "black", fill = "#3B9AB5", shape = 21, size = 4
    ) +
    labs(x = "Mean Intensity", y = "SNR")
}


# ── 9. plotQCOverview ─────────────────────────────────────────────────────────

#' Plot Geary's C vs SNR Overview ('Trashcan Plot')
#'
#' @description 2D scatter of Geary's C vs SNR per marker. A shaded red
#'   rectangle marks the low-quality zone (Geary's C > \code{geary_threshold}
#'   AND SNR < \code{snr_threshold}). Points are labelled using
#'   \pkg{ggrepel}.
#'
#' @param x List output of \code{computeSNR(..., update_correspondence = TRUE)}
#'   where \code{computeGearysC(..., update_correspondence = TRUE)} has also been
#'   run, so \code{CorrespondenceMatrix} contains both \code{GearysC} and
#'   \code{snr} columns.
#' @param geary_threshold Numeric; Geary's C cutoff for the trashcan zone.
#'   Default \code{0.7}.
#' @param snr_threshold Numeric; SNR cutoff for the trashcan zone.
#'   Default \code{3}.
#' @param interactive Logical; if \code{TRUE} returns a \code{plotly} object.
#'   Requires the \pkg{plotly} package. Default \code{FALSE}.
#'
#' @return A \code{ggplot} object, or a \code{plotly} object when
#'   \code{interactive = TRUE}.
#'
#' @examples
#' \dontrun{
#'   processed_qc <- computeGearysC(processed, update_correspondence = TRUE)
#'   processed_qc <- computeSNR(processed_qc, update_correspondence = TRUE)
#'   plotQCOverview(x = processed_qc)
#'   plotQCOverview(x = processed_qc, geary_threshold = 0.8, snr_threshold = 2)
#' }
#' @export
plotQCOverview <- function(x, geary_threshold = 0.7, snr_threshold = 3,
                           interactive = FALSE) {

  if (!requireNamespace("ggrepel", quietly = TRUE)) {
    stop(
      "Package 'ggrepel' is required for plotQCOverview. ",
      "Install it with install.packages('ggrepel')."
    )
  }

  cm <- x$CorrespondenceMatrix
  plot_df <- cm[!is.na(cm$mz_location), ]

  geary_cond  <- plot_df$GearysC > geary_threshold
  snr_cond    <- plot_df$snr < snr_threshold
  plot_df$fill_color <- ifelse(
    geary_cond & snr_cond, "#b51837", "#18B596"
  )

  p <- ggplot() +
    theme_minimal() +
    ggtitle("QC Overview: Geary's C vs SNR") +
    geom_rect(
      aes(xmin = geary_threshold, xmax = 1,
          ymin = 0, ymax = snr_threshold),
      fill = "red3", alpha = 0.5
    ) +
    geom_point(
      data  = plot_df,
      aes(x    = .data[["GearysC"]],
          y    = .data[["snr"]],
          fill = .data[["fill_color"]],
          text = .data[["marker"]]),
      color = "black", shape = 21, size = 4
    ) +
    geom_hline(yintercept = snr_threshold, linetype = "solid", color = "black") +
    geom_vline(xintercept = geary_threshold, linetype = "solid", color = "black") +
    ggrepel::geom_label_repel(
      data          = plot_df,
      aes(x         = .data[["GearysC"]],
          y         = .data[["snr"]],
          label     = .data[["marker"]]),
      color         = "black",
      fill          = "white",
      box.padding   = 0.2,
      point.padding = 0.2,
      max.overlaps  = Inf,
      segment.color = "black",
      size          = 2,
      fontface      = "bold",
      label.size    = 0.2
    ) +
    scale_y_log10() +
    scale_fill_identity() +
    labs(x = "Geary's C", y = "SNR")

  if (!interactive) {
    return(p)
  }

  if (!requireNamespace("plotly", quietly = TRUE)) {
    stop(
      "Package 'plotly' is required for interactive = TRUE. ",
      "Install it with install.packages('plotly')."
    )
  }

  plotly::ggplotly(p, dynamicTicks = FALSE)
}


# ── 10. plotTICHistogram ──────────────────────────────────────────────────────

#' Plot TIC Histogram with Otsu Thresholding
#'
#' @description Histogram of per-pixel total ion count (TIC) on a log10 scale,
#'   with bars coloured by whether pixels fall below (grey = noise) or above
#'   (red = signal) the Otsu threshold.
#'
#' @param x List output of \code{assignMetapeaks}.
#'
#' @return A \code{ggplot} object.
#'
#' @examples
#' \dontrun{
#'   plotTICHistogram(x = processed)
#' }
#' @export
plotTICHistogram <- function(x) {

  tic <- rowSums(x$IntensityDF)
  tic <- tic[tic > 0]

  log_tic   <- log10(tic)
  threshold <- .otsu_thresholding(log_tic)

  plot_df <- data.frame(
    log_tic = log10(1 + tic),
    group   = ifelse(log10(1 + tic) <= threshold, "Noise", "Signal")
  )

  ggplot(plot_df, aes(x = .data[["log_tic"]], fill = .data[["group"]])) +
    geom_histogram(bins = 100, color = "black", alpha = 0.8,
                   position = "identity") +
    geom_vline(xintercept = threshold, linetype = "dashed",
               color = "grey40", linewidth = 1.2) +
    scale_fill_manual(
      values = c("Noise" = "grey", "Signal" = "red3"),
      name   = NULL
    ) +
    theme_minimal() +
    labs(
      title = "Normalised TIC Histogram",
      x     = "log10(TIC)",
      y     = "Count"
    )
}


# ── 11. plotTICSpatial ────────────────────────────────────────────────────────

#' Plot Spatial TIC Signal/Noise Classification
#'
#' @description Spatial scatter plot coloured by TIC-based signal/noise
#'   classification using Otsu thresholding.
#'
#' @param x List output of \code{assignMetapeaks}.
#'
#' @return A \code{ggplot} object.
#'
#' @examples
#' \dontrun{
#'   plotTICSpatial(x = processed)
#' }
#' @export
plotTICSpatial <- function(x) {

  tic      <- rowSums(x$IntensityDF)
  nonzero  <- tic > 0
  tic_nz   <- tic[nonzero]
  coords   <- x$SpatialCoords[nonzero, ]

  log_tic   <- log10(tic_nz)
  threshold <- .otsu_thresholding(log_tic)

  cluster_color <- ifelse(log10(1 + tic_nz) <= threshold, "Noise", "Signal")

  plot_df <- data.frame(coords, cluster_color = cluster_color)

  ggplot(plot_df,
         aes(x = .data[["x"]], y = -.data[["y"]],
             color = factor(.data[["cluster_color"]]))) +
    geom_point(size = 7) +
    scale_color_manual(
      values = c("Noise" = "grey", "Signal" = "red"),
      name   = "Cluster"
    ) +
    theme_void() +
    labs(
      title = "TIC Spatial Distribution",
      x     = "X",
      y     = "Y"
    )
}
