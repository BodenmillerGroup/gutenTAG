test_that("assignMetapeaks works",{

  # get input data
  path <- system.file("extdata/Example_data/Example_data.imzML", package = "gutenTAG")
  panel_path <- system.file("extdata/ref_list.csv", package = "gutenTAG")
  panel <- readPanel(path = panel_path)
  raw <- readMSIData(path)
  pre <- preProcess(raw)
  peaks <- peakDetection(pre)
  metapeaks <- generateMetapeaks(peaks)

  # test that the object after generateMetapeaks is a list
  expect_type(metapeaks, "list")

  # run assignMetapeaks
  cur_test <- assignMetapeaks(metapeaks, pre = pre, refList = panel)

  # test that dimensions are correct
  expect_equal(dim(cur_test$CorrespondenceMatrix), c(13, 7))
  expect_equal(dim(cur_test$IntensityDF), c(64, 13))
  expect_equal(dim(cur_test$SpatialCoords), c(64, 2))
  expect_equal(dim(cur_test$FilteredDF), c(64, 13))
  n_metapeaks <- max(metapeaks$propagation_selection)
  expect_equal(dim(cur_test$AllMetapeaks$AllMetapeaksCorrespondence), c(n_metapeaks, 7))
  expect_equal(dim(cur_test$AllMetapeaks$AllMetapeaksIntensity), c(64, n_metapeaks))
  # check that number of metapeaks matches watershed output
  expect_equal(dim(cur_test$AllMetapeaks$AllMetapeaksIntensity)[2], n_metapeaks)

  # check if output of intensity df is the same
  expect_equal(cur_test$IntensityDF[1][1:10,], c(18.60279186794, 16.0857977316972, 17.9647605811798, 15.0788161493883, 
                                                 16.2613473790588, 9.4926960576303, 14.2638456471353, 10.2109049423084, 
                                                 13.3238938839062, 18.007727014023), tolerance = 0.001)

  # check if spatial coordinates are correct
  expect_equal(range(cur_test$SpatialCoords), c(1, 8))


  # test if character in arguments throws error
  expect_error(assignMetapeaks("test", pre = pre, refList = panel))
  expect_error(assignMetapeaks(metapeaks, pre = "test", refList = panel))
  expect_error(assignMetapeaks(metapeaks, pre = pre, refList = "panel"))

  # test that mz_threshold should be a single numeric
  expect_error(assignMetapeaks(metapeaks, pre = pre, refList = panel, mz_threshold = "test"))

  # test that columns of intensity dataframe are in m/z order
  expect_equal(colnames(cur_test$IntensityDF), cur_test$CorrespondenceMatrix$marker)
  expect_equal(colnames(cur_test$FilteredDF), cur_test$CorrespondenceMatrix$marker)


  # test summarised spectra
  expect_equal(head(cur_test$SummarySpectra$mz),
               c(899.954956054688, 900.010681152344, 900.066345214844, 900.1220703125,
                 900.177795410156, 900.233459472656))
  # skyline spectrum
  expect_equal(head(cur_test$SummarySpectra$skyline),
               c(0, 0.0944778578430998, 0.1707933512843987, 0.4378851361903300,
                 0.3125056625444678, 0.4527385734712326))
  # mean spectrum
  expect_equal(head(cur_test$SummarySpectra$mean),
               c(0, 0.0138080050313242, 0.0190293616937189, 0.0481658760419839,
                 0.0928549807330610, 0.1370896593443408))


})

test_that("params element contains all pipeline parameters with correct default values", {

  path <- system.file("extdata/Example_data/Example_data.imzML", package = "gutenTAG")
  panel_path <- system.file("extdata/ref_list.csv", package = "gutenTAG")
  panel <- readPanel(path = panel_path)
  raw <- readMSIData(path)
  pre <- preProcess(raw)
  peaks <- peakDetection(pre)
  metapeaks <- generateMetapeaks(peaks)
  result <- assignMetapeaks(metapeaks, pre = pre, refList = panel)

  expect_type(result$params, "list")
  expect_named(result$params,
               c("snr", "win", "threshold", "hist_smooth_factor",
                 "fixed.limits", "sparsity", "mz_threshold"))
  expect_equal(result$params$snr, 3)
  expect_equal(result$params$win, 50)
  expect_equal(result$params$threshold, 0.01)
  expect_equal(result$params$hist_smooth_factor, 1)
  expect_null(result$params$fixed.limits)
  expect_equal(result$params$sparsity, 3)
  expect_equal(result$params$mz_threshold, 1)

})

test_that("params captures non-default mz_threshold", {

  path <- system.file("extdata/Example_data/Example_data.imzML", package = "gutenTAG")
  panel_path <- system.file("extdata/ref_list.csv", package = "gutenTAG")
  panel <- readPanel(path = panel_path)
  raw <- readMSIData(path)
  pre <- preProcess(raw)
  peaks <- peakDetection(pre)
  metapeaks <- generateMetapeaks(peaks)
  result <- assignMetapeaks(metapeaks, pre = pre, refList = panel, mz_threshold = 2.5)

  expect_equal(result$params$mz_threshold, 2.5)

})
