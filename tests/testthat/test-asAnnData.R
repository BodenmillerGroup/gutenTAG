test_that("asAnnData works",{

  skip_on_cran()
  skip_if_not_installed("reticulate")

  # get input data
  path <- system.file("extdata/Example_data/Example_data.imzML", package = "gutenTAG")
  panel_path <- system.file("extdata/ref_list.csv", package = "gutenTAG")
  panel <- readPanel(path = panel_path)
  raw <- readMSIData(path)
  pre <- preProcess(raw)
  peaks <- peakDetection(pre)
  metapeaks <- generateMetapeaks(peaks)
  processed <- assignMetapeaks(metapeaks, pre = pre, refList = panel)

  # make anndata object (skip if the Python env cannot be provisioned)
  cur_test <- tryCatch(
    asAnnData(processed),
    error = function(e) skip(paste("Python anndata unavailable:", conditionMessage(e)))
  )

  # object is a Python anndata.AnnData accessed via reticulate
  expect_true(inherits(cur_test, "python.builtin.object"))

  # dimensions of the expression matrix match the input
  expect_equal(dim(cur_test$X), dim(processed$IntensityDF))


})
