test_that("assignMetapeaks works",{

  # get input data
  path <- system.file("extdata/Example_data.imzML", package = "gutenTAG")
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
  expect_equal(dim(cur_test$IntensityDF), c(256, 13))
  expect_equal(dim(cur_test$SpatialCoords), c(256, 2))
  expect_equal(dim(cur_test$FilteredDF), c(256, 13))
  expect_equal(dim(cur_test$AllMetapeaks$AllMetapeaksCorrespondence), c(218, 7))
  expect_equal(dim(cur_test$AllMetapeaks$AllMetapeaksIntensity), c(256, 218))
  # check that number of metapeaks matches watershed output
  expect_equal(dim(cur_test$AllMetapeaks$AllMetapeaksIntensity)[2], max(metapeaks$propagation_selection))

  # check if output of intensity df is the same
  expect_equal(cur_test$IntensityDF[1][1:10,], c(9.7875, 8.8391, 9.2086, 7.8215,
                                                 8.5176, 4.3037, 6.9910, 4.1411,
                                                 5.9847, 4.5863), tolerance = 0.001)

  # check if spatial coordinates are correct
  expect_equal(range(cur_test$SpatialCoords), c(1, 16))


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
               c(0, 0.0944778578430998, 0.189315670352099, 0.43788513619033,
                 0.742802278955215, 0.832960263949027))
  # mean spectrum
  expect_equal(head(cur_test$SummarySpectra$mean),
               c(0, 0.0104144228680699, 0.0239233157355273, 0.0803053741897167,
                 0.155767504458426, 0.203761595637217))


})
