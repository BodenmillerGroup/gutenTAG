test_that("generateMetapeaks works", {
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
  expect_silent(cur_test <- computeVariogram(processed$IntensityDF, coords = processed$SpatialCoords))

  # check if data type is correct
  expect_type(cur_test, "list")

  # check if output is correct
  expect_equal(cur_test[[1]][, 3], c(
    2.72806437546304, 3.13693648542262, 3.42656112710276, 3.4825711368782,
    3.55543933419886, 3.60079212620946, 3.62663470147404, 3.94691460011553,
    3.81095285174737, 4.01796933797204, 3.9819350751035, 3.87224289704557
  ),
  tolerance = 0.001
  )

  # check if output is correct
  expect_equal(cur_test[[1]][, 1], c(
    480, 450, 1288, 1588, 728, 1104, 1010, 1636, 616, 1780, 1088,
    808
  ))


  # test if character in first argument throws error
  expect_error(computeVariogram(df = "processed", coords = processed$SpatialCoords))
  expect_error(computeVariogram(df = processed$IntensityDF, coords = "processed"))
})
