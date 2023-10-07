test_that("computeGMM works",{

  # get input data
  path <- system.file("extdata/Example_data.imzML", package = "maldipackage")
  panel_path <- system.file("extdata/ref_list.csv", package = "maldipackage")
  panel <- readPanel(path = panel_path)
  raw <- readMSIData(path)
  pre <- preProcess(raw, cores = 2)
  peaks <- peakDetection(pre, core = 2)
  metapeaks <- generateMetapeaks(peaks)
  processed <- assignMetapeaks(metapeaks, pre = pre, refList = panel)


  # compute GMM
  cur_test <- computeGMM(x = processed, hist = FALSE)

  # test that output is correct
  expect_equal(cur_test$table[[2]], c(1.4893787954718, 0, 0, 1.42285188639423, 0, 1.69058622560712,
                                      0.470168948733204, 0, 0, 1.570345428177))


  expect_equal(cur_test$table[[3]], c(0.96875, 0, 0, 0.3828125, 0, 0.55078125, 0.078125, 0, 0, 0.609375))

  # test BIC is same
  expect_equal(cur_test$model$BIC[1], 5.4908266)
  # test loglikelihood is the same
  expect_equal(cur_test$model$loglik, 16.608357)


  # test if character in first argument throws error
  expect_error(computeGMM(x = "processed", hist = FALSE))
  computeGMM(x = processed, hist = TRUE)

})
