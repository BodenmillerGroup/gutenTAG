test_that("imageChannel works", {
  # get input data
  path <- system.file("extdata/Example_data.imzML", package = "gutenTAG")
  panel_path <- system.file("extdata/ref_list.csv", package = "gutenTAG")
  panel <- readPanel(path = panel_path)
  raw <- readMSIData(path)
  pre <- preProcess(raw, cores = 2)
  peaks <- peakDetection(pre, core = 2)
  metapeaks <- generateMetapeaks(peaks)
  processed <- assignMetapeaks(metapeaks, pre = pre, refList = panel)


  # compute Geary's C
  expect_silent(cur_test <- imageChannel(x = processed$IntensityDF, coords = processed$SpatialCoords, channel = 1))


  # check if data type is correct
  expect_equal(class(cur_test)[1], "cimg")

  # check if output is correct
  expect_equal(mean(cur_test), 1.891, tolerance = 0.001)

  # check if output is correct
  expect_equal(as.data.frame(cur_test)[, 4][257:280], c(
    6.55226086689995, 4.26364110796174, 2.6607108699066, 6.17412170947713,
    6.20389898040902, 2.93583760073221, 5.2419952122533, 5.30494801453011,
    4.70911068486073, 5.49583488687857, 5.12661053009235, 5.75400695492021,
    3.97242514373091, 10.2807158973856, 5.7218638882349, 7.01733610590797,
    5.60389383788086, 7.9974112204679, 5.07851607822355, 2.91312242777532,
    6.48440505629725, 5.46502750565538, 6.86047612055718, 3.10274556319941
  ),
  tolerance = 0.001
  )


  # test if character in first argument throws error
  expect_error(imageChannel(x = "processed", coords = processed$SpatialCoords, channel_number = 1, quantile_lim = 0.99, interpolate = FALSE, axes = FALSE, colna = FALSE))
  expect_error(imageChannel(x = processed$IntensityDF, coords = "processed", channel_number = 1, quantile_lim = 0.99, interpolate = FALSE, axes = FALSE, colna = FALSE))
  expect_error(imageChannel(x = processed$IntensityDF, coords = processed$SpatialCoords, channel_number = "1", quantile_lim = 0.99, interpolate = FALSE, axes = FALSE, colna = FALSE))
  expect_error(imageChannel(x = processed$IntensityDF, coords = processed$SpatialCoords, channel_number = 1, quantile_lim = "0.99", interpolate = FALSE, axes = FALSE, colna = FALSE))


  ### New tests: additional tests for higher coverage
  # Test with different quantile_lim values
  expect_silent(imageChannel(x = processed$IntensityDF, coords = processed$SpatialCoords, channel_number = 1, quantile_lim = 0))
  expect_silent(imageChannel(x = processed$IntensityDF, coords = processed$SpatialCoords, channel_number = 1, quantile_lim = 1))

  # Test with different channel_number values
  expect_error(imageChannel(x = processed$IntensityDF, coords = processed$SpatialCoords, channel_number = 0))
  expect_error(imageChannel(x = processed$IntensityDF, coords = processed$SpatialCoords, channel_number = ncol(processed$IntensityDF) + 1))

  # Test with different colna values
  expect_silent(imageChannel(x = processed$IntensityDF, coords = processed$SpatialCoords, channel_number = 1, colna = "white"))
  expect_error(imageChannel(x = processed$IntensityDF, coords = processed$SpatialCoords, channel_number = 1, colna = "invalid"))

  # Test with interpolate and axes options
  expect_silent(imageChannel(x = processed$IntensityDF, coords = processed$SpatialCoords, channel_number = 1, interpolate = TRUE))
  expect_silent(imageChannel(x = processed$IntensityDF, coords = processed$SpatialCoords, channel_number = 1, axes = TRUE))
})
