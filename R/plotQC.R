# QC Plot Functions for gutenTAG
# All functions in this file produce quality control plots for MALDI-MSI data.

#' @importFrom ggplot2 ggplot aes annotate geom_bar geom_boxplot geom_histogram
#'   geom_hline geom_line geom_point geom_rect geom_smooth geom_violin geom_vline
#'   scale_color_manual scale_fill_identity scale_fill_manual scale_linetype_manual
#'   scale_x_log10 scale_y_log10 labs ggtitle stat_summary theme theme_minimal
#'   theme_void element_text
#' @importFrom ggbeeswarm geom_quasirandom
#' @importFrom rlang .data
#' @importFrom dplyr arrange filter
#' @importFrom stats lm glm loess sd quantile reorder setNames
#' @importFrom utils modifyList
NULL


# Internal helper for per-layer aesthetic customisation. Builds a ggplot layer
# whose fixed aesthetics come from `defaults`, overridden by anything the user
# supplies in `user` (a named list, e.g. points = list(size = 5, fill = "red")).
# `fixed` holds structural arguments that must NOT be user-overridable
# (data, mapping, yintercept, ...).
.qc_layer <- function(geom, defaults, user = list(), fixed = list()) {
  if (!is.list(user)) {
    stop("Layer aesthetics must be supplied as a named list, ",
         "e.g. points = list(size = 5, fill = \"red\").", call. = FALSE)
  }
  do.call(geom, c(fixed, modifyList(defaults, user)))
}


# Internal helper for recolouring a discrete scale. For layers whose colour/fill
# is mapped to a data grouping (via scale_*_manual), the colours live in the
# scale's palette rather than in a fixed aesthetic. `user` is a named character
# vector overriding entries of `defaults`, e.g. palette = c(Signal = "purple").
.qc_palette <- function(defaults, user = NULL) {
  if (is.null(user)) return(defaults)
  if (!is.character(user) || is.null(names(user))) {
    stop("`palette` must be a named character vector, ",
         "e.g. c(Signal = \"red\", Noise = \"grey\").", call. = FALSE)
  }
  defaults[names(user)] <- user
  defaults
}


# Internal helper: fit the same smoothing model that geom_smooth() draws, so the
# reported statistics match the line. Returns NULL if the method is unknown or
# its package is unavailable (the annotation is then skipped). Mirrors
# geom_smooth()'s "auto" rule (loess for < 1000 points, otherwise gam).
.fit_regression <- function(x, y, method = "lm", formula = y ~ x) {
  df <- data.frame(x = x, y = y)
  if (identical(method, "auto")) {
    method <- if (nrow(df) >= 1000) "gam" else "loess"
  }
  fit_fun <- if (is.function(method)) {
    method
  } else {
    switch(method,
      lm    = lm,
      glm   = glm,
      loess = loess,
      gam   = if (requireNamespace("mgcv", quietly = TRUE)) mgcv::gam,
      rlm   = if (requireNamespace("MASS", quietly = TRUE)) MASS::rlm,
      NULL)
  }
  if (is.null(fit_fun)) return(NULL)
  tryCatch(fit_fun(formula, data = df), error = function(e) NULL)
}


# Internal helper: build the regression stats annotation (a labelled box) for a
# mean-variance-style plot. `fit_x`/`fit_y` are the coordinates in the space the
# line is fitted in (log10 for log scales, raw for linear); `x_pos`/`y_pos` are
# the data-space coordinates to anchor the box. The model is refitted with the
# same method/formula as the drawn line so the reported method, R-squared
# (variance explained) and slope (only when linear in the predictor) describe
# the visible curve. Returns NULL when stats are disabled or cannot be computed.
.regression_stats_layer <- function(fit_x, fit_y, x_pos, y_pos, reg_args,
                                    show_stats, stats) {
  keep <- is.finite(fit_x) & is.finite(fit_y)
  if (!isTRUE(show_stats) || sum(keep) < 2) return(NULL)

  fit <- .fit_regression(fit_x[keep], fit_y[keep], reg_args$method,
                         reg_args$formula)
  if (is.null(fit)) return(NULL)

  yv          <- fit_y[keep]
  fitted_vals <- as.numeric(stats::predict(fit))
  ss_tot      <- sum((yv - mean(yv))^2)
  r2          <- if (ss_tot > 0) 1 - sum((yv - fitted_vals)^2) / ss_tot else NA_real_
  cf          <- tryCatch(stats::coef(fit), error = function(e) NULL)
  slope       <- if (!is.null(cf) && "x" %in% names(cf)) unname(cf[["x"]]) else NA_real_

  method_lbl <- if (is.function(reg_args$method)) {
    "custom"
  } else if (identical(reg_args$method, "auto")) {
    if (sum(keep) >= 1000) "gam" else "loess"
  } else {
    reg_args$method
  }
  parts <- paste0("Method: ", method_lbl)
  if (is.finite(r2))    parts <- c(parts, sprintf("R^2 = %.3f", r2))
  if (is.finite(slope)) parts <- c(parts, sprintf("Slope = %.2f", slope))

  .qc_layer(
    annotate,
    list(x = x_pos, y = y_pos, hjust = 0, vjust = 1, size = 3.5,
         color = "black", fill = "white"),
    stats,
    fixed = list(geom = "label", label = paste(parts, collapse = "\n"))
  )
}


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
#' @param panel Data frame with columns \code{FeatureMass} and \code{Name}.
#'   When loaded via \code{readPanel}, also contains \code{OriginalName}.
#' @param interactive Logical; if \code{TRUE} returns a \code{plotly} object
#'   with mean and skyline spectra on a secondary y-axis. Requires the
#'   \pkg{plotly} package. Default \code{FALSE}.
#'
#' @return A \code{ggplot} object, or a \code{plotly} object when
#'   \code{interactive = TRUE}.
#'
#' @examples
#' rdata_path <- system.file("extdata/Example_data/Example_processed.Rdata", package = "gutenTAG")
#' load(rdata_path)
#' plotMetapeaks(x = results$processed, metapeaks = results$metapeaks,
#'               panel = results$panel)
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

  p <- suppressWarnings(ggplot() +
    theme_minimal() +
    ggtitle("Metapeaks") +
    geom_line(
      data = metapeaks$count_df,
      aes(x = .data[["mz"]], y = .data[["count"]], color = "Count"),
      alpha = 1, linewidth = 0.5
    ) +
    geom_line(
      data = metapeaks$count_smooth_df,
      aes(x = .data[["mz"]], y = .data[["count"]], color = "SmoothCount"),
      alpha = 1, linewidth = 0.5
    ) +
    geom_vline(
      data = panel,
      if (interactive) {
        aes(xintercept = .data[["FeatureMass"]], linetype = "Feature Mass",
            text = .data[["Name"]])
      } else {
        aes(xintercept = .data[["FeatureMass"]], linetype = "Feature Mass")
      },
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
    ))

  # detection threshold (absolute pixel-count) reference line. Guarded so older
  # metapeaks objects without the stored value still plot.
  if (!is.null(metapeaks$params$detection_threshold)) {
    p <- p + geom_hline(
      yintercept = metapeaks$params$detection_threshold,
      linetype = "dashed", color = "grey40", linewidth = 0.6
    )
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
#' @param violin Named list of fixed aesthetics (e.g. \code{fill}, \code{color},
#'   \code{alpha}) for the \code{geom_violin} layer, overriding its defaults.
#' @param boxplot Named list of fixed aesthetics (e.g. \code{fill}, \code{color},
#'   \code{width}) for the \code{geom_boxplot} layer, overriding its defaults.
#'
#' @return A \code{ggplot} object.
#'
#' @examples
#' rdata_path <- system.file("extdata/Example_data/Example_processed.Rdata", package = "gutenTAG")
#' load(rdata_path)
#' plotIntensityDistribution(x = results$processed)
#' plotIntensityDistribution(x = results$processed, violin = list(fill = "skyblue"))
#' @export
plotIntensityDistribution <- function(x, violin = list(), boxplot = list()) {

  df <- x$IntensityDF
  df_long <- utils::stack(df)
  colnames(df_long) <- c("intensity", "marker")

  violin_layer <- .qc_layer(
    geom_violin,
    list(fill = "grey80", color = "black", scale = "width", trim = TRUE),
    violin
  )
  boxplot_layer <- .qc_layer(
    geom_boxplot,
    list(width = 0.15, outlier.shape = NA, fill = "white", color = "black"),
    boxplot
  )

  ggplot(df_long, aes(x = reorder(.data[["marker"]], -.data[["intensity"]]),
                      y = .data[["intensity"]])) +
    violin_layer +
    boxplot_layer +
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
#' @param points Named list of fixed aesthetics (e.g. \code{size}, \code{fill},
#'   \code{color}, \code{alpha}, \code{shape}) for the highlighted/targeted
#'   points layer, overriding its defaults.
#' @param unannotated_points Named list of fixed aesthetics for the grey
#'   unannotated (untargeted) points layer (only drawn when
#'   \code{annotated_only = FALSE}).
#' @param regression Named list of fixed aesthetics (e.g. \code{color},
#'   \code{linewidth}, \code{linetype}, \code{se}) for the regression line,
#'   passed to \code{ggplot2::geom_smooth}. The line is fitted over the
#'   unannotated points, or over the targeted points when
#'   \code{annotated_only = TRUE}. The smoothing method defaults to \code{"lm"};
#'   override it with \code{regression = list(method = ...)}. Available methods
#'   are \code{"lm"}, \code{"glm"}, \code{"gam"} (needs \pkg{mgcv}),
#'   \code{"loess"}, \code{"rlm"} (needs \pkg{MASS}), \code{"auto"} (picks
#'   \code{"loess"} or \code{"gam"} by sample size), or any model function with a
#'   \code{predict} method. The statistics annotation (see \code{show_stats})
#'   refits this same method and reports it, so the box always describes the
#'   visible curve.
#' @param show_stats Logical; if \code{TRUE} (default), annotate the top-left
#'   corner with the fitted model's method, \eqn{R^2} (fraction of variance
#'   explained), and slope (shown only when the model is linear in the
#'   predictor), computed over the same points the line is fitted to.
#' @param stats Named list of fixed aesthetics (e.g. \code{size}, \code{color},
#'   \code{fill}, \code{hjust}, \code{vjust}) for the statistics annotation
#'   label, drawn as a white box with a black outline by default.
#'
#' @return A \code{ggplot} object, or a \code{plotly} object when
#'   \code{interactive = TRUE}.
#'
#' @examples
#' rdata_path <- system.file("extdata/Example_data/Example_processed.Rdata", package = "gutenTAG")
#' load(rdata_path)
#' plotMeanVariance(x = results$processed)
#' plotMeanVariance(x = results$processed, points = list(size = 5, fill = "darkorange"))
#' plotMeanVariance(x = results$processed, regression = list(color = "blue", se = TRUE))
#' plotMeanVariance(x = results$processed, show_stats = FALSE)
#' @export
plotMeanVariance <- function(x, annotated_only = FALSE, interactive = FALSE,
                             points = list(), unannotated_points = list(),
                             regression = list(), show_stats = TRUE,
                             stats = list()) {

  untargeted_df   <- data.frame(x$AllMetapeaks$AllMetapeaksIntensity)
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

  targeted_aes <- if (interactive) {
    aes(x = .data[["targeted_mean"]], y = .data[["targeted_sd"]],
        text = .data[["Name"]])
  } else {
    aes(x = .data[["targeted_mean"]], y = .data[["targeted_sd"]])
  }

  targeted_aes_labeled <- if (interactive) {
    aes(x = .data[["targeted_mean"]], y = .data[["targeted_sd"]],
        text = paste0("Marker: ", .data[["Name"]]))
  } else {
    aes(x = .data[["targeted_mean"]], y = .data[["targeted_sd"]])
  }

  # Build the highlighted/targeted points layer. Per-layer defaults are merged
  # with the user-supplied `points` / `unannotated_points` lists, so any aesthetic
  # can be overridden without enumerating each argument.
  if (annotated_only) {
    targeted_defaults <- list(color = "black", fill = "#3B9AB5", shape = 21, size = 4)
    targeted_mapping  <- targeted_aes
  } else {
    targeted_defaults <- list(color = "black", fill = "red3", shape = 21, size = 3)
    targeted_mapping  <- targeted_aes_labeled
  }
  targeted_layer <- .qc_layer(
    geom_point, targeted_defaults, points,
    fixed = list(data = targeted_stats, mapping = targeted_mapping)
  )

  # The regression line and its stats run over whichever points are plotted:
  # the unannotated (all-metapeak) cloud normally, or the targeted points when
  # annotated_only = TRUE.
  if (annotated_only) {
    reg_stats   <- targeted_stats
    reg_mapping <- aes(x = .data[["targeted_mean"]], y = .data[["targeted_sd"]])
    reg_x       <- targeted_stats$targeted_mean
    reg_y       <- targeted_stats$targeted_sd
  } else {
    reg_stats   <- untargeted_stats
    reg_mapping <- aes(x = .data[["untargeted_mean"]], y = .data[["untargeted_sd"]])
    reg_x       <- untargeted_stats$untargeted_mean
    reg_y       <- untargeted_stats$untargeted_sd
  }

  reg_defaults <- list(method = "lm", formula = y ~ x, color = "red3",
                       se = FALSE, linewidth = 1)
  reg_args <- modifyList(reg_defaults, regression)
  regression_layer <- .qc_layer(
    geom_smooth, reg_defaults, regression,
    fixed = list(data = reg_stats, mapping = reg_mapping)
  )

  # Stats annotation, anchored to the top-left of the plotted points. The fit is
  # in log10 space (matching the log10 scales); passing the log10 coordinates to
  # the helper also drops non-positive values automatically (log10 -> non-finite).
  if (annotated_only) {
    x_all <- targeted_stats$targeted_mean
    y_all <- targeted_stats$targeted_sd
  } else {
    x_all <- c(untargeted_stats$untargeted_mean, targeted_stats$targeted_mean)
    y_all <- c(untargeted_stats$untargeted_sd,   targeted_stats$targeted_sd)
  }
  x_pos <- min(x_all[x_all > 0], na.rm = TRUE)
  y_pos <- max(y_all[y_all > 0], na.rm = TRUE)
  stats_layer <- .regression_stats_layer(
    log10(reg_x), log10(reg_y), x_pos, y_pos, reg_args, show_stats, stats
  )

  if (annotated_only) {
    p <- suppressWarnings(ggplot() +
      theme_minimal() +
      ggtitle("Mean-Variance Plot: Annotated Peaks Only") +
      regression_layer +
      targeted_layer +
      scale_x_log10() +
      scale_y_log10() +
      labs(x = "Mean", y = "SD"))
  } else {
    unannotated_layer <- .qc_layer(
      geom_point,
      list(color = "black", fill = "grey", shape = 21, size = 2),
      unannotated_points,
      fixed = list(
        data    = untargeted_stats,
        mapping = aes(x = .data[["untargeted_mean"]], y = .data[["untargeted_sd"]])
      )
    )
    p <- suppressWarnings(ggplot() +
      theme_minimal() +
      ggtitle("Mean-Variance Plot") +
      unannotated_layer +
      regression_layer +
      targeted_layer +
      scale_x_log10() +
      scale_y_log10() +
      labs(x = "Mean", y = "SD"))
  }

  if (!is.null(stats_layer)) p <- suppressWarnings(p + stats_layer)

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
#' @param points Named list of fixed aesthetics (e.g. \code{size}, \code{fill},
#'   \code{color}, \code{alpha}, \code{shape}) for the highlighted/targeted
#'   points layer, overriding its defaults.
#' @param unannotated_points Named list of fixed aesthetics for the grey
#'   unannotated (untargeted) points layer.
#' @param line Named list of fixed aesthetics (e.g. \code{color},
#'   \code{linetype}, \code{linewidth}) for the zero-reference line.
#'
#' @return A \code{ggplot} object, or a \code{plotly} object when
#'   \code{interactive = TRUE}.
#'
#' @examples
#' rdata_path <- system.file("extdata/Example_data/Example_processed.Rdata", package = "gutenTAG")
#' load(rdata_path)
#' plotMeanVarianceResiduals(x = results$processed)
#' plotMeanVarianceResiduals(x = results$processed,
#'                           points = list(size = 5, fill = "darkorange"),
#'                           line = list(color = "blue", linetype = "dashed"))
#' @export
plotMeanVarianceResiduals <- function(x, standardised = FALSE,
                                      interactive = FALSE, points = list(),
                                      unannotated_points = list(), line = list()) {

  all_df         <- data.frame(x$AllMetapeaks$AllMetapeaksIntensity,
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

  targeted_resid_aes <- if (interactive) {
    aes(x = log(1 + .data[["mean"]]), y = .data[["residual"]],
        text = paste0("Marker: ", .data[["label"]]))
  } else {
    aes(x = log(1 + .data[["mean"]]), y = .data[["residual"]])
  }

  # Per-layer aesthetics: defaults merged with the user-supplied lists.
  targeted_layer <- .qc_layer(
    geom_point, list(color = "black", fill = "red3", shape = 21, size = 3),
    points, fixed = list(data = targeted_df, mapping = targeted_resid_aes)
  )
  unannotated_layer <- .qc_layer(
    geom_point, list(color = "black", fill = "grey", shape = 21, size = 2),
    unannotated_points,
    fixed = list(
      data    = untargeted_df,
      mapping = aes(x = log(1 + .data[["mean"]]), y = .data[["residual"]])
    )
  )
  line_layer <- .qc_layer(
    geom_hline, list(color = "red", linetype = "solid"),
    line, fixed = list(yintercept = 0)
  )

  p <- suppressWarnings(ggplot() +
    theme_minimal() +
    ggtitle(title) +
    line_layer +
    unannotated_layer +
    targeted_layer +
    labs(x = "log(1 + mean)", y = y_label))

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
#' @param bar Named list of fixed aesthetics (e.g. \code{fill}, \code{color},
#'   \code{linewidth}, \code{alpha}) for the \code{geom_bar} layer, overriding
#'   its defaults.
#' @param line Named list of fixed aesthetics (e.g. \code{color},
#'   \code{linetype}, \code{linewidth}) for the mean-score reference line.
#'
#' @return A \code{ggplot} object.
#'
#' @examples
#' rdata_path <- system.file("extdata/Example_data/Example_processed.Rdata", package = "gutenTAG")
#' load(rdata_path)
#' processed_geary <- computeGearysC(results$processed, update_correspondence = TRUE,
#'                                   verbose = FALSE)
#' plotGearysC(x = processed_geary)
#' plotGearysC(x = processed_geary, bar = list(fill = "steelblue"),
#'             line = list(color = "black"))
#' @export
plotGearysC <- function(x, bar = list(), line = list()) {

  cm <- x$CorrespondenceMatrix
  plot_df <- cm[!is.na(cm$GearysC), ]
  plot_df <- plot_df[order(plot_df$GearysC, decreasing = FALSE), ]
  geary_mean <- mean(plot_df$GearysC)

  # Per-layer aesthetics: defaults merged with the user-supplied lists.
  bar_layer <- .qc_layer(
    geom_bar,
    list(stat = "identity", fill = "#b51837", color = "black", linewidth = 0.2),
    bar
  )
  line_layer <- .qc_layer(
    geom_hline, list(linetype = "dashed", color = "grey"),
    line, fixed = list(yintercept = geary_mean)
  )

  ggplot(plot_df, aes(x = reorder(.data[["marker"]], .data[["GearysC"]]),
                      y = .data[["GearysC"]])) +
    bar_layer +
    line_layer +
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
#' @param bar Named list of fixed aesthetics (e.g. \code{color}, \code{width},
#'   \code{alpha}) for the \code{geom_bar} layer. Note that bar \code{fill} is
#'   mapped to the threshold and is controlled via \code{palette}.
#' @param line Named list of fixed aesthetics (e.g. \code{color},
#'   \code{linetype}, \code{linewidth}) for the threshold reference line.
#' @param palette Named character vector overriding the above/below-threshold
#'   fill colours, e.g. \code{c("TRUE" = "firebrick", "FALSE" = "seagreen")}
#'   (\code{"TRUE"} = below threshold).
#'
#' @return A \code{ggplot} object.
#'
#' @examples
#' rdata_path <- system.file("extdata/Example_data/Example_processed.Rdata", package = "gutenTAG")
#' load(rdata_path)
#' processed_snr <- computeSNR(results$processed, update_correspondence = TRUE)
#' plotSNR(x = processed_snr)
#' plotSNR(x = processed_snr, palette = c("TRUE" = "firebrick"))
#' @export
plotSNR <- function(x, snr_threshold = NULL, bar = list(), line = list(),
                    palette = NULL) {

  cm <- x$CorrespondenceMatrix
  plot_df <- cm[!is.na(cm$mz_location), ]
  plot_df <- dplyr::arrange(plot_df, .data[["snr"]])

  if (is.null(snr_threshold)) {
    snr_threshold <- quantile(plot_df$snr, 0.1)
  }

  fill_values <- .qc_palette(c("FALSE" = "#18B596", "TRUE" = "#b51837"), palette)

  bar_layer <- .qc_layer(geom_bar, list(stat = "identity"), bar)
  line_layer <- .qc_layer(
    geom_hline, list(linetype = "dashed", color = "black"),
    line, fixed = list(yintercept = snr_threshold)
  )

  ggplot(plot_df,
         aes(x = reorder(.data[["marker"]], -.data[["snr"]]),
             y = .data[["snr"]],
             fill = .data[["snr"]] <= snr_threshold)) +
    bar_layer +
    scale_fill_manual(
      values = fill_values,
      labels = c("FALSE" = "Above threshold", "TRUE" = "Below threshold"),
      name   = NULL
    ) +
    line_layer +
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
#' @param histogram Named list of fixed aesthetics (e.g. \code{bins},
#'   \code{color}, \code{alpha}) for the \code{geom_histogram} layer. Note that
#'   fill is mapped to the signal/noise cluster and is set via \code{palette}.
#' @param palette Named character vector overriding the cluster fill colours,
#'   e.g. \code{c(Signal = "purple", Noise = "grey80")}.
#'
#' @return A \code{ggplot} object.
#'
#' @examples
#' rdata_path <- system.file("extdata/Example_data/Example_processed.Rdata", package = "gutenTAG")
#' load(rdata_path)
#' processed_snr <- computeSNR(results$processed, update_correspondence = TRUE)
#' plotSNRHistogram(x = processed_snr, channel = 1)
#' plotSNRHistogram(x = processed_snr, channel = 1, palette = c(Signal = "purple"))
#' @export
plotSNRHistogram <- function(x, channel, histogram = list(), palette = NULL) {

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

  fill_values <- .qc_palette(c("Noise" = "grey", "Signal" = "red"), palette)

  hist_layer <- .qc_layer(
    geom_histogram,
    list(bins = 100, color = "black", alpha = 0.5, position = "identity"),
    histogram
  )

  ggplot(hist_df, aes(x = .data[["count"]], fill = .data[["cluster_color"]])) +
    hist_layer +
    scale_fill_manual(
      values = fill_values,
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
#' @param points Named list of fixed aesthetics (e.g. \code{size}, \code{shape},
#'   \code{alpha}) for the \code{geom_point} layer. Note that point colour is
#'   mapped to the signal/noise cluster and is set via \code{palette}.
#' @param palette Named character vector overriding the cluster colours,
#'   e.g. \code{c(Signal = "purple", Noise = "grey80")}.
#'
#' @return A \code{ggplot} object.
#'
#' @examples
#' rdata_path <- system.file("extdata/Example_data/Example_processed.Rdata", package = "gutenTAG")
#' load(rdata_path)
#' processed_snr <- computeSNR(results$processed, update_correspondence = TRUE)
#' plotSpatialSNR(x = processed_snr, channel = 1)
#' plotSpatialSNR(x = processed_snr, channel = 1, points = list(size = 4))
#' @export
plotSpatialSNR <- function(x, channel, points = list(), palette = NULL) {

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

  color_values <- .qc_palette(c("Noise" = "grey", "Signal" = "red"), palette)
  point_layer <- .qc_layer(geom_point, list(size = 7), points)

  ggplot(plot_df,
         aes(x = .data[["x"]], y = -.data[["y"]],
             color = factor(.data[["cluster_color"]]))) +
    point_layer +
    scale_color_manual(
      values = color_values,
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
#' @param points Named list of fixed aesthetics (e.g. \code{size}, \code{fill},
#'   \code{color}, \code{alpha}, \code{shape}) for the \code{geom_point} layer,
#'   overriding its defaults.
#' @param regression Named list of fixed aesthetics (e.g. \code{color},
#'   \code{linewidth}, \code{linetype}, \code{se}) for the regression line,
#'   passed to \code{ggplot2::geom_smooth}. The smoothing method defaults to
#'   \code{"lm"}; override it with \code{regression = list(method = ...)}.
#'   Available methods are \code{"lm"}, \code{"glm"}, \code{"gam"} (needs
#'   \pkg{mgcv}), \code{"loess"}, \code{"rlm"} (needs \pkg{MASS}), \code{"auto"},
#'   or any model function with a \code{predict} method. The statistics
#'   annotation (see \code{show_stats}) refits this same method and reports it.
#' @param show_stats Logical; if \code{TRUE} (default), annotate the top-left
#'   corner with the fitted model's method, \eqn{R^2} (fraction of variance
#'   explained), and slope (shown only when the model is linear in the
#'   predictor).
#' @param stats Named list of fixed aesthetics (e.g. \code{size}, \code{color},
#'   \code{fill}, \code{hjust}, \code{vjust}) for the statistics annotation
#'   label, drawn as a white box with a black outline by default.
#'
#' @return A \code{ggplot} object.
#'
#' @examples
#' rdata_path <- system.file("extdata/Example_data/Example_processed.Rdata", package = "gutenTAG")
#' load(rdata_path)
#' processed_snr <- computeSNR(results$processed, update_correspondence = TRUE)
#' plotMeanVsSNR(x = processed_snr)
#' plotMeanVsSNR(x = processed_snr, points = list(size = 6, fill = "darkorange"))
#' plotMeanVsSNR(x = processed_snr, regression = list(method = "loess"))
#' @export
plotMeanVsSNR <- function(x, points = list(), regression = list(),
                          show_stats = TRUE, stats = list()) {

  targeted_mean <- colMeans(x$IntensityDF)
  targeted_snr  <- x$CorrespondenceMatrix$snr
  plot_df <- data.frame(
    targeted_mean = targeted_mean,
    targeted_snr  = targeted_snr,
    Name          = colnames(x$IntensityDF)
  )

  snr_aes <- aes(x = .data[["targeted_mean"]], y = .data[["targeted_snr"]])

  # Per-layer aesthetics: defaults merged with the user-supplied `points` list.
  point_layer <- .qc_layer(
    geom_point, list(color = "black", fill = "#3B9AB5", shape = 21, size = 4),
    points, fixed = list(data = plot_df, mapping = snr_aes)
  )

  reg_defaults <- list(method = "lm", formula = y ~ x, color = "red3",
                       se = FALSE, linewidth = 1)
  reg_args <- modifyList(reg_defaults, regression)
  regression_layer <- .qc_layer(
    geom_smooth, reg_defaults, regression,
    fixed = list(data = plot_df, mapping = snr_aes)
  )

  # Axes are linear here, so the stats model is fitted on the raw mean/SNR
  # values (no log transform) to match the drawn line.
  finite <- is.finite(targeted_mean) & is.finite(targeted_snr)
  x_pos  <- if (any(finite)) min(targeted_mean[finite]) else NA_real_
  y_pos  <- if (any(finite)) max(targeted_snr[finite])  else NA_real_
  stats_layer <- .regression_stats_layer(
    targeted_mean, targeted_snr, x_pos, y_pos, reg_args, show_stats, stats
  )

  p <- ggplot() +
    theme_minimal() +
    ggtitle("Mean-SNR Plot") +
    point_layer +
    regression_layer +
    labs(x = "Mean Intensity", y = "SNR")

  if (!is.null(stats_layer)) p <- p + stats_layer
  p
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
#' @param points Named list of fixed aesthetics (e.g. \code{size}, \code{shape})
#'   for the \code{geom_point} layer. Point fill encodes the quality zone and is
#'   set via \code{palette}.
#' @param rect Named list of fixed aesthetics (e.g. \code{fill}, \code{alpha})
#'   for the shaded low-quality zone rectangle.
#' @param hline,vline Named lists of fixed aesthetics (e.g. \code{color},
#'   \code{linetype}) for the SNR (horizontal) and Geary's C (vertical)
#'   threshold lines.
#' @param labels Named list of fixed aesthetics passed to
#'   \code{ggrepel::geom_label_repel} for the marker labels.
#' @param palette Named character vector overriding the point zone colours,
#'   with names \code{"low_quality"} and \code{"pass"}, e.g.
#'   \code{c(low_quality = "firebrick")}.
#'
#' @return A \code{ggplot} object, or a \code{plotly} object when
#'   \code{interactive = TRUE}.
#'
#' @examples
#' rdata_path <- system.file("extdata/Example_data/Example_processed.Rdata", package = "gutenTAG")
#' load(rdata_path)
#' processed_qc <- computeGearysC(results$processed, update_correspondence = TRUE)
#' processed_qc <- computeSNR(processed_qc, update_correspondence = TRUE)
#' plotQCOverview(x = processed_qc)
#' plotQCOverview(x = processed_qc, points = list(size = 5))
#' @export
plotQCOverview <- function(x, geary_threshold = 0.7, snr_threshold = 3,
                           interactive = FALSE, points = list(), rect = list(),
                           hline = list(), vline = list(), labels = list(),
                           palette = NULL) {

  if (!requireNamespace("ggrepel", quietly = TRUE)) {
    stop(
      "Package 'ggrepel' is required for plotQCOverview. ",
      "Install it with install.packages('ggrepel')."
    )
  }

  cm <- x$CorrespondenceMatrix
  plot_df <- cm[!is.na(cm$mz_location), ]

  zone_colors <- .qc_palette(
    c(low_quality = "#b51837", pass = "#18B596"), palette
  )
  geary_cond  <- plot_df$GearysC > geary_threshold
  snr_cond    <- plot_df$snr < snr_threshold
  plot_df$fill_color <- ifelse(
    geary_cond & snr_cond, zone_colors[["low_quality"]], zone_colors[["pass"]]
  )

  point_mapping <- if (interactive) {
    aes(x    = .data[["GearysC"]], y = .data[["snr"]],
        fill = .data[["fill_color"]], text = .data[["marker"]])
  } else {
    aes(x    = .data[["GearysC"]], y = .data[["snr"]],
        fill = .data[["fill_color"]])
  }

  rect_layer <- .qc_layer(
    geom_rect, list(fill = "red3", alpha = 0.5), rect,
    fixed = list(mapping = aes(xmin = geary_threshold, xmax = 1,
                               ymin = 0, ymax = snr_threshold))
  )
  point_layer <- .qc_layer(
    geom_point, list(color = "black", shape = 21, size = 4), points,
    fixed = list(data = plot_df, mapping = point_mapping)
  )
  hline_layer <- .qc_layer(
    geom_hline, list(linetype = "solid", color = "black"), hline,
    fixed = list(yintercept = snr_threshold)
  )
  vline_layer <- .qc_layer(
    geom_vline, list(linetype = "solid", color = "black"), vline,
    fixed = list(xintercept = geary_threshold)
  )
  label_layer <- .qc_layer(
    ggrepel::geom_label_repel,
    list(color = "black", fill = "white", box.padding = 0.2,
         point.padding = 0.2, max.overlaps = Inf, segment.color = "black",
         size = 2, fontface = "bold", label.size = 0.2),
    labels,
    fixed = list(
      data    = plot_df,
      mapping = aes(x = .data[["GearysC"]], y = .data[["snr"]],
                    label = .data[["marker"]])
    )
  )

  p <- suppressWarnings(ggplot() +
    theme_minimal() +
    ggtitle("QC Overview: Geary's C vs SNR") +
    rect_layer +
    point_layer +
    hline_layer +
    vline_layer +
    label_layer +
    scale_y_log10() +
    scale_fill_identity() +
    labs(x = "Geary's C", y = "SNR"))

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
#' @param histogram Named list of fixed aesthetics (e.g. \code{bins},
#'   \code{color}, \code{alpha}) for the \code{geom_histogram} layer. Fill is
#'   mapped to the signal/noise group and is set via \code{palette}.
#' @param vline Named list of fixed aesthetics (e.g. \code{color},
#'   \code{linetype}, \code{linewidth}) for the Otsu-threshold line.
#' @param palette Named character vector overriding the group fill colours,
#'   e.g. \code{c(Signal = "navy", Noise = "grey80")}.
#'
#' @return A \code{ggplot} object.
#'
#' @examples
#' rdata_path <- system.file("extdata/Example_data/Example_processed.Rdata", package = "gutenTAG")
#' load(rdata_path)
#' plotTICHistogram(x = results$processed)
#' plotTICHistogram(x = results$processed, palette = c(Signal = "navy"))
#' @export
plotTICHistogram <- function(x, histogram = list(), vline = list(),
                             palette = NULL) {

  tic <- rowSums(x$IntensityDF)
  tic <- tic[tic > 0]

  log_tic   <- log10(tic)
  threshold <- .otsu_thresholding(log_tic)

  plot_df <- data.frame(
    log_tic = log10(1 + tic),
    group   = ifelse(log10(1 + tic) <= threshold, "Noise", "Signal")
  )

  fill_values <- .qc_palette(c("Noise" = "grey", "Signal" = "red3"), palette)

  hist_layer <- .qc_layer(
    geom_histogram,
    list(bins = 100, color = "black", alpha = 0.8, position = "identity"),
    histogram
  )
  vline_layer <- .qc_layer(
    geom_vline, list(linetype = "dashed", color = "grey40", linewidth = 1.2),
    vline, fixed = list(xintercept = threshold)
  )

  ggplot(plot_df, aes(x = .data[["log_tic"]], fill = .data[["group"]])) +
    hist_layer +
    vline_layer +
    scale_fill_manual(
      values = fill_values,
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
#' @param points Named list of fixed aesthetics (e.g. \code{size}, \code{shape},
#'   \code{alpha}) for the \code{geom_point} layer. Point colour is mapped to the
#'   signal/noise group and is set via \code{palette}.
#' @param palette Named character vector overriding the group colours,
#'   e.g. \code{c(Signal = "navy", Noise = "grey80")}.
#'
#' @return A \code{ggplot} object.
#'
#' @examples
#' rdata_path <- system.file("extdata/Example_data/Example_processed.Rdata", package = "gutenTAG")
#' load(rdata_path)
#' plotTICSpatial(x = results$processed)
#' plotTICSpatial(x = results$processed, points = list(size = 4))
#' @export
plotTICSpatial <- function(x, points = list(), palette = NULL) {

  tic      <- rowSums(x$IntensityDF)
  nonzero  <- tic > 0
  tic_nz   <- tic[nonzero]
  coords   <- x$SpatialCoords[nonzero, ]

  log_tic   <- log10(tic_nz)
  threshold <- .otsu_thresholding(log_tic)

  cluster_color <- ifelse(log10(1 + tic_nz) <= threshold, "Noise", "Signal")

  plot_df <- data.frame(coords, cluster_color = cluster_color)

  color_values <- .qc_palette(c("Noise" = "grey", "Signal" = "red"), palette)
  point_layer <- .qc_layer(geom_point, list(size = 7), points)

  ggplot(plot_df,
         aes(x = .data[["x"]], y = -.data[["y"]],
             color = factor(.data[["cluster_color"]]))) +
    point_layer +
    scale_color_manual(
      values = color_values,
      name   = "Cluster"
    ) +
    theme_void() +
    labs(
      title = "TIC Spatial Distribution",
      x     = "X",
      y     = "Y"
    )
}


# ── 12. plotMassShift ─────────────────────────────────────────────────────────

#' Plot Mass Shift Between Expected and Picked m/z
#'
#' @description Visualises the deviation between each targeted peak's expected
#'   m/z (the panel feature mass) and its picked m/z (the metapeak max), computed
#'   as \code{mz_location - expected_mz_location}. Two views are available:
#'   \code{"distribution"} draws a single-column quasirandom (beeswarm) plot of
#'   the mass-shift distribution across all targeted peaks, with the mean marked
#'   by a diamond; \code{"spectrum"} plots the same mass shift against the
#'   expected m/z, so any systematic drift across the mass range is visible. A
#'   grey reference line is drawn at zero shift.
#'
#' @param x List output of \code{\link{assignMetapeaks}}. Its
#'   \code{CorrespondenceMatrix} must contain \code{mz_location} and
#'   \code{expected_mz_location} columns.
#' @param type Character; \code{"distribution"} (default) for the beeswarm of the
#'   mass-shift distribution, or \code{"spectrum"} for mass shift along the m/z
#'   axis.
#' @param points Named list of fixed aesthetics (e.g. \code{size}, \code{fill},
#'   \code{color}, \code{alpha}, \code{shape}) for the points layer (the
#'   quasirandom points for \code{"distribution"}, the scatter for
#'   \code{"spectrum"}).
#' @param mean_point Named list of fixed aesthetics for the mean-marker diamond
#'   (\code{"distribution"} view only).
#' @param line Named list of fixed aesthetics (e.g. \code{color},
#'   \code{linewidth}, \code{linetype}) for the zero-reference line.
#'
#' @return A \code{ggplot} object.
#'
#' @examples
#' rdata_path <- system.file("extdata/Example_data/Example_processed.Rdata", package = "gutenTAG")
#' load(rdata_path)
#' plotMassShift(x = results$processed, type = "distribution")
#' plotMassShift(x = results$processed, type = "spectrum")
#' plotMassShift(x = results$processed, points = list(fill = "darkorange"))
#' @export
plotMassShift <- function(x, type = c("distribution", "spectrum"),
                          points = list(), mean_point = list(), line = list()) {

  type <- match.arg(type)

  cm <- x$CorrespondenceMatrix
  df <- cm[!is.na(cm$mz_location) & !is.na(cm$expected_mz_location), ]
  df$diff <- df$mz_location - df$expected_mz_location   # Picked - Expected (m/z)

  line_layer <- .qc_layer(
    geom_hline, list(color = "grey60", linewidth = 0.6),
    line, fixed = list(yintercept = 0)
  )

  if (type == "distribution") {

    df$group <- "All peaks"

    points_layer <- .qc_layer(
      geom_quasirandom,
      list(shape = 21, stroke = 0.2, width = 0.25, alpha = 0.4, size = 4,
           color = "black", fill = "#3B9AB5"),
      points
    )
    mean_layer <- .qc_layer(
      stat_summary,
      list(shape = 23, size = 4, fill = "black", color = "black", stroke = 1),
      mean_point, fixed = list(fun = mean, geom = "point")
    )

    ggplot(df, aes(x = .data[["group"]], y = .data[["diff"]])) +
      points_layer +
      mean_layer +
      line_layer +
      theme_minimal() +
      labs(
        x     = NULL,
        y     = "Picked - Expected (m/z)",
        title = "Mass Shift Distribution"
      )

  } else {  # "spectrum": shift along the m/z axis

    points_layer <- .qc_layer(
      geom_point,
      list(shape = 21, size = 3, color = "black", fill = "#3B9AB5", alpha = 0.7),
      points
    )

    ggplot(df, aes(x = .data[["expected_mz_location"]], y = .data[["diff"]])) +
      points_layer +
      line_layer +
      theme_minimal() +
      labs(
        x     = "Expected m/z",
        y     = "Picked - Expected (m/z)",
        title = "Mass Shift Across m/z Range"
      )

  }
}
