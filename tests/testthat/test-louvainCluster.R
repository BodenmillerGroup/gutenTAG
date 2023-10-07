test_that("louvainCluster works",{

  # get input data
  path <- system.file("extdata/Example_data.imzML", package = "maldipackage")
  panel_path <- system.file("extdata/ref_list.csv", package = "maldipackage")
  panel <- readPanel(path = panel_path)
  raw <- readMSIData(path)
  pre <- preProcess(raw, cores = 2)
  peaks <- peakDetection(pre, core = 2)
  metapeaks <- generateMetapeaks(peaks)
  processed <- assignMetapeaks(metapeaks, pre = pre, refList = panel)

  # extract data
  expect_type(df <- processed$IntensityDF, "list")
  expect_type(coords <- processed$SpatialCoords, "list")

  # compute louvainCluster
  expect_silent(cur_test <- louvainCluster(x = df, coords = coords, k = 30 , metric = "angular", resolution = 0.5))
  expect_warning(louvainCluster(x = df, coords = coords))

  # test that output is correct is correct
  expect_equal(length(cur_test), 256)
  expect_equal(as.numeric(cur_test[1:10]), c(1, 1, 1, 2, 2, 3, 3, 3, 3, 2))


  # test if character in first argument throws error
  expect_error(louvainCluster(x = "df", coords = coords, k = 30 , metric = "angular", resolution = 0.5))
  expect_error(louvainCluster(x = df, coords = "coords", k = 30 , metric = "angular", resolution = 0.5))
  expect_error(louvainCluster(x = df, coords = coords, k = "30" , metric = "angular", resolution = 0.5))
  expect_error(louvainCluster(x = df, coords = coords, k = 30 , metric = "angular", resolution = "0.5"))

  # expect warning
  expect_warning(louvainCluster(processed$IntensityDF, processed$SpatialCoords))


})
