test_that("generateMetapeaks works",{

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
  expect_equal(cur_test[[1]][, 3], c(2.08189234849682, 2.40373103483664, 2.72038977319223, 2.75868424897782,
                                     2.78745223346144, 2.80389531814698, 2.86367048369323, 3.13216729674623,
                                     2.98872350399834, 3.18468042493146, 3.16670828148136, 3.04492277141782))

  # check if output is correct
  expect_equal(cur_test[[1]][, 1], c(480, 450, 1288, 1588, 728, 1104, 1010, 1636, 616, 1780, 1088,
                                     808))


  # test if character in first argument throws error
  expect_error(computeVariogram(df = "processed", coords = processed$SpatialCoords))
  expect_error(computeVariogram(df = processed$IntensityDF, coords = "processed"))


})
