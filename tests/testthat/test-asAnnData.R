test_that("asAnnData works",{

  # get input data
  path <- system.file("extdata/Example_data.imzML", package = "gutenTAG")
  panel_path <- system.file("extdata/ref_list.csv", package = "gutenTAG")
  panel <- readPanel(path = panel_path)
  raw <- readMSIData(path)
  pre <- preProcess(raw, cores = 2)
  peaks <- peakDetection(pre, core = 2)
  metapeaks <- generateMetapeaks(peaks)
  processed <- assignMetapeaks(metapeaks, pre = pre, refList = panel)

  # make anndata object
  cur_test <- asAnnData(processed)

  # test that object inherits from R6 class and AnnData class
  expect_true(inherits(cur_test, "R6"))
  expect_true(inherits(cur_test, "AnnDataR6"))

  # test that dimensions of expression matrix in anndata object are the same as the input
  expect_equal(dim(cur_test$X), dim(processed$IntensityDF))


})
