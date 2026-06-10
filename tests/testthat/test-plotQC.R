test_that("plotMetapeaks works", {

  p <- plotMetapeaks(x = qc_processed, metapeaks = qc_metapeaks, panel = qc_panel)
  expect_s3_class(p, "ggplot")
  expect_equal(p$labels$title, "Metapeaks")

  skip_if_not_installed("plotly")
  p_interactive <- plotMetapeaks(x = qc_processed, metapeaks = qc_metapeaks,
                                 panel = qc_panel, interactive = TRUE)
  expect_s3_class(p_interactive, "plotly")

})

test_that("plotMetapeaks draws the detection-threshold line when available", {

  has_hline <- function(p) {
    any(vapply(p$layers, function(l) inherits(l$geom, "GeomHline"), logical(1)))
  }

  # inject the stored value (bundled fixture predates it) -> line present
  mp <- qc_metapeaks
  mp$params$detection_threshold <- 5
  expect_true(has_hline(plotMetapeaks(x = qc_processed, metapeaks = mp, panel = qc_panel)))

  # absent value -> guarded, no line drawn (backward compatible)
  mp_no <- qc_metapeaks
  mp_no$params$detection_threshold <- NULL
  expect_false(has_hline(plotMetapeaks(x = qc_processed, metapeaks = mp_no, panel = qc_panel)))

})

test_that("plotMassShift works", {

  p_dist <- plotMassShift(x = qc_processed, type = "distribution")
  expect_s3_class(p_dist, "ggplot")
  expect_equal(p_dist$labels$title, "Mass Shift Distribution")
  expect_equal(p_dist$labels$y, "Picked - Expected (m/z)")

  p_spec <- plotMassShift(x = qc_processed, type = "spectrum")
  expect_s3_class(p_spec, "ggplot")
  expect_equal(p_spec$labels$title, "Mass Shift Across m/z Range")
  expect_equal(p_spec$labels$x, "Expected m/z")

  # default type is distribution
  expect_equal(plotMassShift(x = qc_processed)$labels$title, "Mass Shift Distribution")

  # invalid type errors via match.arg
  expect_error(plotMassShift(x = qc_processed, type = "nonsense"))

  # per-layer lists: points (layer 1) and zero line (layer 3)
  p_custom <- plotMassShift(x = qc_processed, points = list(fill = "darkorange"),
                            line = list(color = "blue"))
  expect_equal(p_custom$layers[[1]]$aes_params$fill, "darkorange")
  expect_equal(p_custom$layers[[3]]$aes_params$colour, "blue")

})

test_that("plotIntensityDistribution works", {

  p <- plotIntensityDistribution(x = qc_processed)
  expect_s3_class(p, "ggplot")
  expect_equal(p$labels$title, "Mean Marker Intensities")
  expect_equal(p$labels$y, "Intensity per pixel")

  # violin (layer 1) and boxplot (layer 2) per-layer lists
  p_custom <- plotIntensityDistribution(x = qc_processed,
                                        violin = list(fill = "skyblue"),
                                        boxplot = list(fill = "wheat"))
  expect_equal(p_custom$layers[[1]]$aes_params$fill, "skyblue")
  expect_equal(p_custom$layers[[2]]$aes_params$fill, "wheat")

})

test_that("plotMeanVariance works", {

  p <- plotMeanVariance(x = qc_processed)
  expect_s3_class(p, "ggplot")
  expect_equal(p$labels$title, "Mean-Variance Plot")

  p_annotated <- plotMeanVariance(x = qc_processed, annotated_only = TRUE)
  expect_s3_class(p_annotated, "ggplot")
  expect_equal(p_annotated$labels$title, "Mean-Variance Plot: Annotated Peaks Only")

  # Default-mode layer order: 1 unannotated, 2 regression, 3 targeted, 4 stats
  # `points` list sets fixed aesthetics on the highlighted/targeted layer (3)
  p_custom <- plotMeanVariance(x = qc_processed,
                               points = list(size = 5, fill = "darkorange"))
  expect_equal(p_custom$layers[[3]]$aes_params$size, 5)
  expect_equal(p_custom$layers[[3]]$aes_params$fill, "darkorange")

  # `unannotated_points` list sets fixed aesthetics on the grey layer (1)
  p_bg <- plotMeanVariance(x = qc_processed, unannotated_points = list(fill = "khaki"))
  expect_equal(p_bg$layers[[1]]$aes_params$fill, "khaki")

  # regression line (layer 2) drawn by default and customisable
  expect_s3_class(p$layers[[2]]$geom, "GeomSmooth")
  expect_equal(p$layers[[2]]$aes_params$colour, "red3")
  p_reg <- plotMeanVariance(x = qc_processed, regression = list(color = "blue"))
  expect_equal(p_reg$layers[[2]]$aes_params$colour, "blue")

  # annotated_only also draws a regression + stats, fitted over the targeted
  # points: layers are regression (1), targeted (2), stats (3)
  p_ann <- plotMeanVariance(x = qc_processed, annotated_only = TRUE)
  expect_length(p_ann$layers, 3)
  expect_s3_class(p_ann$layers[[1]]$geom, "GeomSmooth")
  expect_s3_class(p_ann$layers[[3]]$geom, "GeomLabel")
  # show_stats = FALSE drops the annotation in annotated_only mode too
  expect_length(plotMeanVariance(x = qc_processed, annotated_only = TRUE,
                                 show_stats = FALSE)$layers, 2)

  # targeted defaults preserved when nothing is passed (layer 3)
  expect_equal(p$layers[[3]]$aes_params$fill, "red3")
  expect_equal(p$layers[[3]]$aes_params$size, 3)

  # stats annotation (layer 4) is a labelled box, toggleable and styleable
  expect_s3_class(p$layers[[4]]$geom, "GeomLabel")
  expect_equal(p$layers[[4]]$aes_params$fill, "white")
  expect_equal(p$layers[[4]]$aes_params$colour, "black")
  expect_length(plotMeanVariance(x = qc_processed, show_stats = FALSE)$layers, 3)
  p_stats <- plotMeanVariance(x = qc_processed, stats = list(size = 6, fill = "lightyellow"))
  expect_equal(p_stats$layers[[4]]$aes_params$size, 6)
  expect_equal(p_stats$layers[[4]]$aes_params$fill, "lightyellow")

  # the annotation must resolve to finite coords so it renders under log scales
  # (a -Inf anchor would become NaN and silently vanish)
  txt <- suppressWarnings(ggplot2::layer_data(p, 4))
  expect_true(all(is.finite(c(txt$x, txt$y))))
  txt_ann <- suppressWarnings(
    ggplot2::layer_data(plotMeanVariance(x = qc_processed, annotated_only = TRUE), 3))
  expect_true(all(is.finite(c(txt_ann$x, txt_ann$y))))

  # stats follow the regression method and label it; slope only for linear fits
  expect_match(p$layers[[4]]$aes_params$label, "Method: lm")
  expect_match(p$layers[[4]]$aes_params$label, "Slope")
  p_loess <- plotMeanVariance(x = qc_processed, regression = list(method = "loess"))
  expect_match(p_loess$layers[[4]]$aes_params$label, "Method: loess")
  expect_no_match(p_loess$layers[[4]]$aes_params$label, "Slope")

  skip_if_not_installed("plotly")
  p_interactive <- plotMeanVariance(x = qc_processed, interactive = TRUE)
  expect_s3_class(p_interactive, "plotly")

})

test_that("plotMeanVarianceResiduals works", {

  p <- plotMeanVarianceResiduals(x = qc_processed)
  expect_s3_class(p, "ggplot")
  expect_equal(p$labels$title, "Residuals vs Mean Intensity")
  expect_equal(p$labels$y, "Residuals")

  p_std <- plotMeanVarianceResiduals(x = qc_processed, standardised = TRUE)
  expect_s3_class(p_std, "ggplot")
  expect_equal(p_std$labels$title, "Standardised Residuals vs Mean Intensity")
  expect_equal(p_std$labels$y, "Standardised residuals")

  # `points` list sets fixed aesthetics on the highlighted/targeted (last) layer
  p_custom <- plotMeanVarianceResiduals(x = qc_processed,
                                        points = list(size = 5, fill = "darkorange"))
  resid_layer <- p_custom$layers[[length(p_custom$layers)]]
  expect_equal(resid_layer$aes_params$size, 5)
  expect_equal(resid_layer$aes_params$fill, "darkorange")

  # `line` list controls the zero-reference line (first layer); default is red
  expect_equal(p$layers[[1]]$aes_params$colour, "red")
  p_line <- plotMeanVarianceResiduals(x = qc_processed,
                                      line = list(color = "blue", linetype = "dashed"))
  expect_equal(p_line$layers[[1]]$aes_params$colour, "blue")
  expect_equal(p_line$layers[[1]]$aes_params$linetype, "dashed")

  skip_if_not_installed("plotly")
  p_interactive <- plotMeanVarianceResiduals(x = qc_processed, interactive = TRUE)
  expect_s3_class(p_interactive, "plotly")

})

test_that("plotGearysC works", {

  p <- plotGearysC(x = qc_processed_geary)
  expect_s3_class(p, "ggplot")
  expect_equal(p$labels$title, "Geary's C Score")
  expect_equal(p$labels$y, "Geary's C")

  # `bar` / `line` lists set fixed aesthetics on their respective layers
  p_custom <- plotGearysC(x = qc_processed_geary,
                          bar = list(fill = "steelblue"),
                          line = list(color = "black"))
  expect_equal(p_custom$layers[[1]]$aes_params$fill, "steelblue")
  expect_equal(p_custom$layers[[2]]$aes_params$colour, "black")

  # all-NA GearysC: filter removes all rows, ggplot should still build
  x_all_na <- qc_processed_geary
  x_all_na$CorrespondenceMatrix$GearysC <- NA_real_
  expect_s3_class(plotGearysC(x = x_all_na), "ggplot")

})

test_that("plotSNR works", {

  p <- plotSNR(x = qc_processed_snr)
  expect_s3_class(p, "ggplot")
  expect_equal(p$labels$title, "SNR per Marker")
  expect_equal(p$labels$y, "SNR")

  p_thresh <- plotSNR(x = qc_processed_snr, snr_threshold = 5)
  expect_s3_class(p_thresh, "ggplot")

  # `bar` list sets non-mapped aesthetics; `palette` recolours; both build
  p_custom <- plotSNR(x = qc_processed_snr, bar = list(alpha = 0.5),
                      palette = c("TRUE" = "firebrick"))
  expect_s3_class(p_custom, "ggplot")
  expect_equal(p_custom$layers[[1]]$aes_params$alpha, 0.5)

})

test_that("plotSNRHistogram works", {

  chan_name <- colnames(qc_processed_snr$IntensityDF)[1]

  p_int  <- plotSNRHistogram(x = qc_processed_snr, channel = 1)
  p_char <- plotSNRHistogram(x = qc_processed_snr, channel = chan_name)
  expect_s3_class(p_int, "ggplot")
  expect_s3_class(p_char, "ggplot")

  # integer and character resolve to the same channel
  expect_equal(p_int$labels$title, p_char$labels$title)
  expect_equal(p_int$labels$title, paste("Histogram of Signal and Noise:", chan_name))

  # `histogram` list + `palette` recolour build and apply
  p_custom <- plotSNRHistogram(x = qc_processed_snr, channel = 1,
                               histogram = list(alpha = 0.3),
                               palette = c(Signal = "purple"))
  expect_s3_class(p_custom, "ggplot")
  expect_equal(p_custom$layers[[1]]$aes_params$alpha, 0.3)

})

test_that("plotSpatialSNR works", {

  chan_name <- colnames(qc_processed_snr$IntensityDF)[1]

  p_int  <- plotSpatialSNR(x = qc_processed_snr, channel = 1)
  p_char <- plotSpatialSNR(x = qc_processed_snr, channel = chan_name)
  expect_s3_class(p_int, "ggplot")
  expect_s3_class(p_char, "ggplot")
  expect_equal(p_int$labels$title, paste("Spatial Distribution:", chan_name))

  # `points` list + `palette` recolour build and apply
  p_custom <- plotSpatialSNR(x = qc_processed_snr, channel = 1,
                             points = list(size = 4),
                             palette = c(Signal = "purple"))
  expect_s3_class(p_custom, "ggplot")
  expect_equal(p_custom$layers[[1]]$aes_params$size, 4)

})

test_that("plotMeanVsSNR works", {

  p <- plotMeanVsSNR(x = qc_processed_snr)
  expect_s3_class(p, "ggplot")
  expect_equal(p$labels$title, "Mean-SNR Plot")
  expect_equal(p$labels$x, "Mean Intensity")
  expect_equal(p$labels$y, "SNR")

  # `points` list sets fixed aesthetics on the points layer (layer 1)
  p_custom <- plotMeanVsSNR(x = qc_processed_snr,
                            points = list(size = 6, fill = "darkorange"))
  pts_layer <- p_custom$layers[[1]]
  expect_equal(pts_layer$aes_params$size, 6)
  expect_equal(pts_layer$aes_params$fill, "darkorange")

  # regression line (layer 2) + method-coupled stats box (layer 3)
  expect_s3_class(p$layers[[2]]$geom, "GeomSmooth")
  expect_s3_class(p$layers[[3]]$geom, "GeomLabel")
  expect_equal(p$layers[[3]]$aes_params$fill, "white")
  expect_match(p$layers[[3]]$aes_params$label, "Method: lm")
  expect_match(p$layers[[3]]$aes_params$label, "Slope")
  txt <- suppressWarnings(ggplot2::layer_data(p, 3))
  expect_true(all(is.finite(c(txt$x, txt$y))))

  # method couples to the stats box; loess omits slope
  p_loess <- plotMeanVsSNR(x = qc_processed_snr, regression = list(method = "loess"))
  expect_match(p_loess$layers[[3]]$aes_params$label, "Method: loess")
  expect_no_match(p_loess$layers[[3]]$aes_params$label, "Slope")

  # show_stats = FALSE drops the annotation (points + regression only)
  expect_length(plotMeanVsSNR(x = qc_processed_snr, show_stats = FALSE)$layers, 2)

})

test_that("plotQCOverview works", {

  p <- plotQCOverview(x = qc_processed_snr)
  expect_s3_class(p, "ggplot")
  expect_equal(p$labels$x, "Geary's C")
  expect_equal(p$labels$y, "SNR")

  p_custom <- plotQCOverview(x = qc_processed_snr,
                             geary_threshold = 0.8, snr_threshold = 2)
  expect_s3_class(p_custom, "ggplot")

  # `points` list (layer 2, after the rect) + `palette` recolour build and apply
  p_aes <- plotQCOverview(x = qc_processed_snr, points = list(size = 5),
                          palette = c(low_quality = "firebrick"))
  expect_s3_class(p_aes, "ggplot")
  expect_equal(p_aes$layers[[2]]$aes_params$size, 5)

  skip_if_not_installed("plotly")
  p_interactive <- plotQCOverview(x = qc_processed_snr, interactive = TRUE)
  expect_s3_class(p_interactive, "plotly")

})

test_that("plotTICHistogram works", {

  p <- plotTICHistogram(x = qc_processed)
  expect_s3_class(p, "ggplot")
  expect_equal(p$labels$title, "Normalised TIC Histogram")
  expect_equal(p$labels$x, "log10(TIC)")
  expect_equal(p$labels$y, "Count")

  # histogram (layer 1) + vline (layer 2) lists, plus `palette`
  p_custom <- plotTICHistogram(x = qc_processed, histogram = list(alpha = 0.3),
                               vline = list(color = "black"),
                               palette = c(Signal = "navy"))
  expect_s3_class(p_custom, "ggplot")
  expect_equal(p_custom$layers[[1]]$aes_params$alpha, 0.3)
  expect_equal(p_custom$layers[[2]]$aes_params$colour, "black")

})

test_that("plotTICSpatial works", {

  p <- plotTICSpatial(x = qc_processed)
  expect_s3_class(p, "ggplot")
  expect_equal(p$labels$title, "TIC Spatial Distribution")

  n_nonzero <- sum(rowSums(qc_processed$IntensityDF) > 0)
  expect_equal(nrow(ggplot2::layer_data(p)), n_nonzero)

  # `points` list + `palette` recolour build and apply
  p_custom <- plotTICSpatial(x = qc_processed, points = list(size = 4),
                             palette = c(Signal = "navy"))
  expect_s3_class(p_custom, "ggplot")
  expect_equal(p_custom$layers[[1]]$aes_params$size, 4)

})

