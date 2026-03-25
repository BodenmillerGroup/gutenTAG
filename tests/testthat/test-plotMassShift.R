# tests/testthat/test-plotMassShift.R
#
# Tests for computePixelMassShift() and the four plotMassShift*() functions.
#
# Shared fixtures (qc_processed, qc_metapeaks, qc_panel) are provided by
# helper-plotQC.R and are available automatically in every test block.
#
# The 'pre' object (MSImagingExperiment from preProcess()) is not included in
# Example_processed.Rdata, so it is constructed once here at file scope using
# the pattern established in test-computeSNR.R.

.imzml_path <- system.file("extdata/Example_data.imzML", package = "gutenTAG")
.raw_pre    <- Cardinal::readMSIData(.imzml_path)
pre         <- gutenTAG::preProcess(.raw_pre)

# Full pipeline result used by computePixelMassShift tests and plot tests that
# require $PixelMassShift.
qc_shift <- computePixelMassShift(
  qc_processed,
  pre                   = pre,
  metapeaks             = qc_metapeaks,
  update_correspondence = TRUE
)


# ── 1. computePixelMassShift ──────────────────────────────────────────────────

test_that("computePixelMassShift returns correct structure when update_correspondence = FALSE", {

  result <- computePixelMassShift(
    qc_processed,
    pre       = pre,
    metapeaks = qc_metapeaks
  )

  expect_true(is.list(result))
  expect_named(result, c("MassShiftDF", "PixelMassShift"), ignore.order = FALSE)

  # MassShiftDF: one row per pixel, one column per observed marker
  n_pixels   <- nrow(qc_processed$SpatialCoords)
  n_observed <- sum(!is.na(qc_processed$CorrespondenceMatrix$mz_location))
  expect_true(is.data.frame(result$MassShiftDF))
  expect_equal(nrow(result$MassShiftDF), n_pixels)
  expect_equal(ncol(result$MassShiftDF), n_observed)

  # PixelMassShift: numeric vector of length n_pixels; all-zero pixels become NA
  expect_true(is.numeric(result$PixelMassShift))
  expect_equal(length(result$PixelMassShift), n_pixels)

  # All-zero pixels must be NA, not NaN
  na_pixels <- is.na(result$PixelMassShift)
  expect_false(any(is.nan(result$PixelMassShift[na_pixels])))
})

test_that("computePixelMassShift returns augmented list when update_correspondence = TRUE", {

  result <- computePixelMassShift(
    qc_processed,
    pre                   = pre,
    metapeaks             = qc_metapeaks,
    update_correspondence = TRUE
  )

  # Must return the full x list, not just MassShiftDF/PixelMassShift
  expect_true(is.list(result))
  expect_true("CorrespondenceMatrix" %in% names(result))
  expect_true("IntensityDF"          %in% names(result))
  expect_true("MassShiftDF"          %in% names(result))
  expect_true("PixelMassShift"       %in% names(result))

  # mean_mass_shift column added to CorrespondenceMatrix
  expect_true("mean_mass_shift" %in% colnames(result$CorrespondenceMatrix))

  # Unobserved markers (NA mz_location) get NA mean_mass_shift
  unobs_rows <- is.na(result$CorrespondenceMatrix$mz_location)
  expect_true(all(is.na(result$CorrespondenceMatrix$mean_mass_shift[unobs_rows])))

  # Observed markers must not have NaN mean_mass_shift
  obs_shifts <- result$CorrespondenceMatrix$mean_mass_shift[!unobs_rows]
  expect_false(any(is.nan(obs_shifts)))
})

test_that("computePixelMassShift errors on invalid x", {

  expect_error(
    computePixelMassShift("not_a_list", pre = pre, metapeaks = qc_metapeaks),
    regexp = "'x' should be a list"
  )

  x_missing_field <- qc_processed[setdiff(names(qc_processed), "SpatialCoords")]
  expect_error(
    computePixelMassShift(x_missing_field, pre = pre, metapeaks = qc_metapeaks),
    regexp = "'x' should be a list"
  )
})

test_that("computePixelMassShift errors when pre is not an MSImagingExperiment", {

  expect_error(
    computePixelMassShift(qc_processed, pre = "not_pre", metapeaks = qc_metapeaks),
    regexp = "'pre' must be an MSImagingExperiment"
  )
})

test_that("computePixelMassShift errors when metapeaks is missing required fields", {

  bad_metapeaks <- list(metapeaks = list(max = qc_metapeaks$metapeaks$max))
  expect_error(
    computePixelMassShift(qc_processed, pre = pre, metapeaks = bad_metapeaks),
    regexp = "'metapeaks' must be a list"
  )
})

test_that("computePixelMassShift errors when all mz_location are NA", {

  x_all_na <- qc_processed
  x_all_na$CorrespondenceMatrix$mz_location <- NA_real_

  expect_error(
    computePixelMassShift(x_all_na, pre = pre, metapeaks = qc_metapeaks),
    regexp = "No observed markers found"
  )
})

test_that("computePixelMassShift errors when update_correspondence is not logical", {

  expect_error(
    computePixelMassShift(
      qc_processed, pre = pre, metapeaks = qc_metapeaks,
      update_correspondence = "yes"
    ),
    regexp = "'update_correspondence' should be a boolean"
  )
})


# ── 2. plotMassShiftDistribution ──────────────────────────────────────────────

test_that("plotMassShiftDistribution returns a ggplot with the correct title", {

  p <- plotMassShiftDistribution(x = qc_processed)

  expect_s3_class(p, "ggplot")
  expect_equal(p$labels$title, "Mass Shift Distribution")
})

test_that("plotMassShiftDistribution label_threshold controls ggrepel layer", {

  # Very large threshold: no points qualify, ggrepel layer absent
  p_no_labels <- plotMassShiftDistribution(x = qc_processed, label_threshold = 1e9)
  layer_classes_no <- vapply(p_no_labels$layers, function(l) class(l$geom)[1], character(1))
  expect_false(any(grepl("GeomTextRepel", layer_classes_no)))

  # Zero threshold: every point qualifies, ggrepel layer present
  p_all_labels <- plotMassShiftDistribution(x = qc_processed, label_threshold = 0)
  layer_classes_all <- vapply(p_all_labels$layers, function(l) class(l$geom)[1], character(1))
  expect_true(any(grepl("GeomTextRepel", layer_classes_all)))
})

test_that("plotMassShiftDistribution interactive = TRUE returns plotly", {

  skip_if_not_installed("plotly")
  p <- plotMassShiftDistribution(x = qc_processed, interactive = TRUE)
  expect_s3_class(p, "plotly")
})

test_that("plotMassShiftDistribution errors on invalid x", {

  expect_error(plotMassShiftDistribution(x = "not_a_list"))
  expect_error(plotMassShiftDistribution(x = list()))
})


# ── 3. plotMassShiftVsMz ──────────────────────────────────────────────────────

test_that("plotMassShiftVsMz returns a ggplot with the correct title and axis labels", {

  p <- plotMassShiftVsMz(x = qc_processed)

  expect_s3_class(p, "ggplot")
  expect_equal(p$labels$title, "Mass Shift vs m/z")
  expect_equal(p$labels$x, "Expected m/z (Da)")
  expect_equal(p$labels$y, "Mass shift (Da)")
})

test_that("plotMassShiftVsMz interactive = TRUE returns plotly", {

  skip_if_not_installed("plotly")
  p <- plotMassShiftVsMz(x = qc_processed, interactive = TRUE)
  expect_s3_class(p, "plotly")
})

test_that("plotMassShiftVsMz errors on invalid x", {

  expect_error(plotMassShiftVsMz(x = "not_a_list"))
  expect_error(plotMassShiftVsMz(x = list()))
})


# ── 4. plotMassShiftVsIntensity ───────────────────────────────────────────────

test_that("plotMassShiftVsIntensity returns a ggplot with the correct title and axis labels", {

  p <- plotMassShiftVsIntensity(x = qc_processed)

  expect_s3_class(p, "ggplot")
  expect_equal(p$labels$title, "Mass Shift vs Mean Intensity")
  expect_equal(p$labels$x, "log10(Mean Intensity + 1)")
  expect_equal(p$labels$y, "Mass shift (Da)")
})

test_that("plotMassShiftVsIntensity interactive = TRUE returns plotly", {

  skip_if_not_installed("plotly")
  p <- plotMassShiftVsIntensity(x = qc_processed, interactive = TRUE)
  expect_s3_class(p, "plotly")
})

test_that("plotMassShiftVsIntensity errors on invalid x", {

  expect_error(plotMassShiftVsIntensity(x = "not_a_list"))
  expect_error(plotMassShiftVsIntensity(x = list()))
})


# ── 5. plotMassShiftSpatial ───────────────────────────────────────────────────

test_that("plotMassShiftSpatial returns a ggplot with the correct title", {

  p <- plotMassShiftSpatial(x = qc_shift)

  expect_s3_class(p, "ggplot")
  expect_equal(p$labels$title, "Spatial Mass Shift")
})

test_that("plotMassShiftSpatial plots only non-NA pixels", {

  n_nonzero <- sum(!is.na(qc_shift$PixelMassShift))
  p <- plotMassShiftSpatial(x = qc_shift)
  expect_equal(nrow(ggplot2::layer_data(p)), n_nonzero)
})

test_that("plotMassShiftSpatial errors when PixelMassShift is absent", {

  x_no_shift <- qc_shift
  x_no_shift$PixelMassShift <- NULL

  expect_error(
    plotMassShiftSpatial(x = x_no_shift),
    regexp = "computePixelMassShift"
  )
})

test_that("plotMassShiftSpatial errors when all PixelMassShift values are NA", {

  x_all_na <- qc_shift
  x_all_na$PixelMassShift <- rep(NA_real_, length(qc_shift$PixelMassShift))

  expect_error(
    plotMassShiftSpatial(x = x_all_na),
    regexp = "No pixels with a computed mass shift"
  )
})

test_that("plotMassShiftSpatial works from a synthetic minimal x", {

  # Confirm that the function requires only $PixelMassShift and $SpatialCoords
  x_minimal <- list(
    PixelMassShift = c(-0.02, 0.00, 0.01, -0.01, 0.03, NA, 0.02, -0.03, 0.00),
    SpatialCoords  = data.frame(x = rep(1:3, each = 3), y = rep(1:3, times = 3))
  )

  p <- plotMassShiftSpatial(x = x_minimal)

  expect_s3_class(p, "ggplot")
  # Eight non-NA pixels should be plotted
  expect_equal(nrow(ggplot2::layer_data(p)), 8L)
})
