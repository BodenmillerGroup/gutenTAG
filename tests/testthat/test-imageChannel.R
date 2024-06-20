test_that("imageChannel works",{

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
  expect_equal(mean(cur_test), 1.59034635)

  # check if outut is correct
  expect_equal(as.data.frame(cur_test)[, 4][257:300], c(4.76116810842076, 3.17607185993165, 1.92611508943016, 4.56589048246013,
                                                        5.20620718577363, 2.70269539185386, 3.74961860055951, 4.4526247100637,
                                                        3.79729292427373, 4.13266727212503, 4.47581109283329, 4.88015570193393,
                                                        3.04574901373604, 8.44431899054, 4.76621901143297, 6.30595739216694,
                                                        4.39780775471946, 6.02808179747598, 3.7412976139573, 2.30679575419345,
                                                        4.97865551762262, 2.79422234909605, 5.49652130707497, 2.63360638530468,
                                                        5.33746713005442, 4.35070529983597, 3.75030975371312, 4.52395325680042,
                                                        4.10985903783833, 4.56964765207574, 4.73855038299939, 1.53321451804249,
                                                        4.25333444385425, 3.08503620779209, 3.22723680947978, 5.41591989324583,
                                                        2.57613747769858, 3.037671376565, 4.36592172155993, 4.65419595571034,
                                                        4.67545491562303, 5.16613269676679, 4.46922536268614, 2.67442868003569))


  # test if character in first argument throws error
  expect_error(imageChannel(x = "processed", coords = processed$SpatialCoords, channel_number = 1, quantile_lim = 0.99, interpolate = FALSE, axes = FALSE, colna = FALSE))
  expect_error(imageChannel(x = processed$IntensityDF, coords = "processed", channel_number = 1, quantile_lim = 0.99, interpolate = FALSE, axes = FALSE, colna = FALSE))
  expect_error(imageChannel(x = processed$IntensityDF, coords = processed$SpatialCoords, channel_number = "1", quantile_lim = 0.99, interpolate = FALSE, axes = FALSE, colna = FALSE))
  expect_error(imageChannel(x = processed$IntensityDF, coords = processed$SpatialCoords, channel_number = 1, quantile_lim = "0.99", interpolate = FALSE, axes = FALSE, colna = FALSE))


})
