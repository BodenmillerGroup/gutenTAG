test_that("computeSNR works",{

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
  cur_test <- computeSNR(processed, update_correspondence = TRUE)

  # test if character in first argument throws error
  expect_error(computeSNR("processed"))

  # Test correspondence matrix
  expect_true("snr" %in% colnames(cur_test$CorrespondenceMatrix))
  expect_equal(cur_test$CorrespondenceMatrix$snr, c(1, 1.46662215920315, 0, 0, 2.00872721014486, 0, 3.08641724395242,
                                                    3.80602499772131, 3.01605464957345, 15.2477686105868, 2.38579865682754,
                                                    1.41332759446791, 2.9132994473942))

})

