# Tests for the `transformation` argument of plotIntensityDistribution.
# The default "none" behaviour is covered in test-plotQC.R; these tests focus on
# the log / z-scaled transformations merged in from the former
# plotMarkerMeanIntensity.

test_that("plotIntensityDistribution returns a ggplot for each transformation", {

  for (tr in c("none", "log", "z-scaled")) {
    p <- plotIntensityDistribution(x = qc_processed, transformation = tr)
    expect_s3_class(p, "ggplot")
    expect_s3_class(p$layers[[1]]$geom, "GeomViolin")
    expect_s3_class(p$layers[[2]]$geom, "GeomBoxplot")
  }
})

test_that("plotIntensityDistribution y-label reflects the transformation", {

  expect_equal(
    plotIntensityDistribution(x = qc_processed, transformation = "none")$labels$y,
    "Intensity per pixel"
  )
  expect_equal(
    plotIntensityDistribution(x = qc_processed, transformation = "log")$labels$y,
    "Intensity per pixel (log10)"
  )
  expect_equal(
    plotIntensityDistribution(x = qc_processed,
                              transformation = "z-scaled")$labels$y,
    "Intensity per pixel (log10 z-score)"
  )
})

test_that("plotIntensityDistribution errors on an invalid transformation", {

  expect_error(
    plotIntensityDistribution(x = qc_processed, transformation = "nonsense")
  )
})

test_that("plotIntensityDistribution z-scaled values are centred and scaled", {

  df_long <- utils::stack(qc_processed$IntensityDF)
  colnames(df_long) <- c("intensity", "marker")
  expected <- as.numeric(scale(log10(df_long$intensity + 1)))

  p   <- plotIntensityDistribution(x = qc_processed, transformation = "z-scaled")
  vals <- p$data$value

  # pooled across all pixels/markers: global mean 0, sd 1
  expect_equal(mean(vals), 0, tolerance = 1e-6)
  expect_equal(sd(vals),   1, tolerance = 1e-6)
  expect_equal(sort(vals), sort(expected), tolerance = 1e-6)
})

test_that("plotIntensityDistribution log transform is finite with zero-intensity pixels", {

  # the fixture's per-pixel intensities include zeros; log10(x + 1) must stay finite
  expect_true(any(qc_processed$IntensityDF == 0))

  p <- plotIntensityDistribution(x = qc_processed, transformation = "log")
  expect_true(all(is.finite(p$data$value)))
})

test_that("plotIntensityDistribution marker x-order is identical across transformations", {

  # reorder(marker, -intensity) sets the discrete x order; read it back from the
  # built plot's x-axis labels (layer_data maps x to integer positions only)
  levels_of <- function(tr) {
    p  <- plotIntensityDistribution(x = qc_processed, transformation = tr)
    pb <- ggplot2::ggplot_build(p)
    pb$layout$panel_params[[1]]$x$get_labels()
  }

  none_lv <- levels_of("none")
  expect_identical(none_lv, levels_of("log"))
  expect_identical(none_lv, levels_of("z-scaled"))
})
