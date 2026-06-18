test_that("plotMassShiftVsIntensity works", {

  p <- plotMassShiftVsIntensity(x = qc_processed)
  expect_s3_class(p, "ggplot")
  expect_equal(p$labels$title, "Mass Shift vs Intensity")
  expect_equal(p$labels$x, "Picked - Expected (m/z)")
  expect_equal(p$labels$y, "Mean intensity (log10)")

  # layer 1 is the zero-reference line (vertical, at shift = 0), layer 2 the points
  expect_s3_class(p$layers[[1]]$geom, "GeomVline")
  expect_s3_class(p$layers[[2]]$geom, "GeomPoint")

  # y axis is log10: the rendered y positions are log10 of the mean intensity
  cm  <- qc_processed$CorrespondenceMatrix
  df  <- cm[!is.na(cm$mz_location) & !is.na(cm$expected_mz_location) &
              !is.na(cm$mean) & cm$mean > 0, ]
  pts <- ggplot2::layer_data(p, 2)
  expect_equal(sort(pts$y), sort(log10(df$mean)), tolerance = 1e-8)

  # slightly portrait aspect ratio (>= 1)
  expect_equal(p$theme$aspect.ratio, 1.2)
})

test_that("plotMassShiftVsIntensity keeps one point per qualifying marker", {

  cm <- qc_processed$CorrespondenceMatrix
  n_expected <- sum(!is.na(cm$mz_location) & !is.na(cm$expected_mz_location) &
                      !is.na(cm$mean) & cm$mean > 0)

  p   <- plotMassShiftVsIntensity(x = qc_processed)
  pts <- ggplot2::layer_data(p, 2)
  expect_equal(nrow(pts), n_expected)
})

test_that("plotMassShiftVsIntensity drops NA / non-positive intensities with a warning", {

  x_bad <- qc_processed
  cm    <- x_bad$CorrespondenceMatrix
  # force a couple of markers to be undefined on a log axis
  cm$mean[1] <- NA_real_
  cm$mean[2] <- 0
  x_bad$CorrespondenceMatrix <- cm

  expect_warning(plotMassShiftVsIntensity(x = x_bad), regexp = "non-positive")

  n_expected <- sum(!is.na(cm$mz_location) & !is.na(cm$expected_mz_location) &
                      !is.na(cm$mean) & cm$mean > 0)
  p   <- suppressWarnings(plotMassShiftVsIntensity(x = x_bad))
  pts <- ggplot2::layer_data(p, 2)
  expect_equal(nrow(pts), n_expected)
})

test_that("plotMassShiftVsIntensity colours points by shift direction", {

  cm <- qc_processed$CorrespondenceMatrix
  df <- cm[!is.na(cm$mz_location) & !is.na(cm$expected_mz_location) &
             !is.na(cm$mean) & cm$mean > 0, ]
  diff <- df$mz_location - df$expected_mz_location

  p     <- plotMassShiftVsIntensity(x = qc_processed, show_labels = FALSE)
  fills <- ggplot2::layer_data(p, 2)$fill

  # below the zero line (negative shift) is blue, above (>= 0) is red.
  # checked by count to stay independent of layer row order
  expect_equal(sum(fills == "#b51837"), sum(diff >= 0))   # above -> red
  expect_equal(sum(fills == "#3B9AB5"), sum(diff <  0))   # below -> blue
})

test_that("plotMassShiftVsIntensity respects palette override", {

  p <- plotMassShiftVsIntensity(x = qc_processed, show_labels = FALSE,
                                palette = c(above = "#000000", below = "#ffffff"))
  fills <- ggplot2::layer_data(p, 2)$fill
  expect_true(all(fills %in% c("#000000", "#ffffff")))
})

test_that("plotMassShiftVsIntensity ignores `fill` in points so palette still applies", {

  # a stray fill in `points` must not clobber the direction colouring
  expect_warning(
    p <- plotMassShiftVsIntensity(x = qc_processed, show_labels = FALSE,
                                  points = list(size = 6, fill = "orange"),
                                  palette = c(above = "red", below = "blue")),
    regexp = "fill.*ignored"
  )
  fills <- ggplot2::layer_data(p, 2)$fill
  expect_false("orange" %in% fills)            # the stray fill is dropped
  expect_true(all(fills %in% c("red", "blue"))) # palette colours win
})

test_that("plotMassShiftVsIntensity errors on invalid input", {

  expect_error(plotMassShiftVsIntensity(x = list()),
               regexp = "CorrespondenceMatrix")
})

test_that("plotMassShiftVsIntensity honours user point aesthetics", {

  p <- plotMassShiftVsIntensity(x = qc_processed, points = list(size = 6),
                                show_labels = FALSE)
  expect_equal(p$layers[[2]]$aes_params$size, 6)
})
