test_that("plotMetapeaks works", {

  p <- plotMetapeaks(x = qc_processed, metapeaks = qc_metapeaks, panel = qc_panel)
  expect_s3_class(p, "ggplot")
  expect_equal(p$labels$title, "Metapeaks")

  skip_if_not_installed("plotly")
  p_interactive <- plotMetapeaks(x = qc_processed, metapeaks = qc_metapeaks,
                                 panel = qc_panel, interactive = TRUE)
  expect_s3_class(p_interactive, "plotly")

})

test_that("plotIntensityDistribution works", {

  p <- plotIntensityDistribution(x = qc_processed)
  expect_s3_class(p, "ggplot")
  expect_equal(p$labels$title, "Mean Marker Intensities")
  expect_equal(p$labels$y, "Intensity per pixel")

})

test_that("plotMeanVariance works", {

  p <- plotMeanVariance(x = qc_processed)
  expect_s3_class(p, "ggplot")
  expect_equal(p$labels$title, "Mean-Variance Plot")

  p_annotated <- plotMeanVariance(x = qc_processed, annotated_only = TRUE)
  expect_s3_class(p_annotated, "ggplot")
  expect_equal(p_annotated$labels$title, "Mean-Variance Plot: Annotated Peaks Only")

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

  skip_if_not_installed("plotly")
  p_interactive <- plotMeanVarianceResiduals(x = qc_processed, interactive = TRUE)
  expect_s3_class(p_interactive, "plotly")

})

test_that("plotGearysC works", {

  p <- plotGearysC(x = qc_processed_geary)
  expect_s3_class(p, "ggplot")
  expect_equal(p$labels$title, "Geary's C Score")
  expect_equal(p$labels$y, "Geary's C")

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

})

test_that("plotSpatialSNR works", {

  chan_name <- colnames(qc_processed_snr$IntensityDF)[1]

  p_int  <- plotSpatialSNR(x = qc_processed_snr, channel = 1)
  p_char <- plotSpatialSNR(x = qc_processed_snr, channel = chan_name)
  expect_s3_class(p_int, "ggplot")
  expect_s3_class(p_char, "ggplot")
  expect_equal(p_int$labels$title, paste("Spatial Distribution:", chan_name))

})

test_that("plotMeanVsSNR works", {

  p <- plotMeanVsSNR(x = qc_processed_snr)
  expect_s3_class(p, "ggplot")
  expect_equal(p$labels$title, "Mean-SNR Plot")
  expect_equal(p$labels$x, "Mean Intensity")
  expect_equal(p$labels$y, "SNR")

})

test_that("plotQCOverview works", {

  p <- plotQCOverview(x = qc_processed_snr)
  expect_s3_class(p, "ggplot")
  expect_equal(p$labels$x, "Geary's C")
  expect_equal(p$labels$y, "SNR")

  p_custom <- plotQCOverview(x = qc_processed_snr,
                             geary_threshold = 0.8, snr_threshold = 2)
  expect_s3_class(p_custom, "ggplot")

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

})

test_that("plotTICSpatial works", {

  p <- plotTICSpatial(x = qc_processed)
  expect_s3_class(p, "ggplot")
  expect_equal(p$labels$title, "TIC Spatial Distribution")

  n_nonzero <- sum(rowSums(qc_processed$IntensityDF) > 0)
  expect_equal(nrow(ggplot2::layer_data(p)), n_nonzero)

})

