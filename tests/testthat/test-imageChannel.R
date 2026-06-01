test_that("imageChannel works", {

  # get input data
  path <- system.file("extdata/Example_data/Example_data.imzML", package = "gutenTAG")
  panel_path <- system.file("extdata/ref_list.csv", package = "gutenTAG")
  panel <- readPanel(path = panel_path)
  raw <- readMSIData(path)
  pre <- preProcess(raw)
  peaks <- peakDetection(pre)
  metapeaks <- generateMetapeaks(peaks)
  processed <- assignMetapeaks(metapeaks, pre = pre, refList = panel)


  # test that imageChannel returns a ggplot object
  expect_s3_class(imageChannel(processed, channel = 1), "gg")


  # test list (assignMetapeaks) input — by index and by name
  expect_silent(imageChannel(processed, channel = 1))
  expect_silent(imageChannel(processed, channel = colnames(processed$IntensityDF)[1]))


  # test plain data.frame + coords input
  expect_silent(imageChannel(processed$IntensityDF, coords = processed$SpatialCoords, channel = 1))


  # test quantile thresholding — max fill value equals max raw intensity when threshold = 1.0
  intensity <- processed$IntensityDF[, 1]
  p_full <- imageChannel(processed, channel = 1, quantile_threshold = 1.0)
  expect_equal(max(p_full$data$intensity, na.rm = TRUE), max(intensity, na.rm = TRUE))

  # test quantile thresholding — max fill value is capped at median when threshold = 0.5
  p_half <- imageChannel(processed, channel = 1, quantile_threshold = 0.5)
  expect_equal(max(p_half$data$intensity, na.rm = TRUE), median(intensity, na.rm = TRUE))


  # test palette argument
  expect_silent(imageChannel(processed, channel = 1, palette = "magma"))
  expect_error(imageChannel(processed, channel = 1, palette = 123))


  # test na_colour argument
  expect_silent(imageChannel(processed, channel = 1, na_colour = "white"))
  expect_error(imageChannel(processed, channel = 1, na_colour = TRUE))


  # test invalid channel
  expect_error(imageChannel(processed, channel = "nonexistent_marker"), "not found")
  expect_error(imageChannel(processed, channel = TRUE))


  # test invalid quantile_threshold
  expect_error(imageChannel(processed, channel = 1, quantile_threshold = 0))
  expect_error(imageChannel(processed, channel = 1, quantile_threshold = 1.5))
  expect_error(imageChannel(processed, channel = 1, quantile_threshold = "high"))


  # test that coords is required when x is a data.frame
  expect_error(imageChannel(processed$IntensityDF, coords = NULL))


  # test that theme arguments are passed through via ...
  expect_silent(imageChannel(processed, channel = 1, plot.title = ggplot2::element_text(size = 20)))


})
