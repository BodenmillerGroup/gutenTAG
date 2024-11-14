test_that("computeSNR works", {
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
  cur_test <- computeSNR(processed, update_correspondence = TRUE, package = "mclust")

  # test if character in first argument throws error
  expect_error(computeSNR("processed"))

  # Test correspondence matrix
  expect_true("snr" %in% colnames(cur_test$CorrespondenceMatrix))
  expect_equal(cur_test$CorrespondenceMatrix$snr, c(
    2.07919158184831, 1.48489368089413, 0, 0, 2.22915785723627,
    0, 3.1144754240626, 3.72397173486482, 3.323669305973, 16.9176896774671,
    2.40921027994661, 1.41160848191838, 3.03517059385205
  ))

  # test flexmix snr
  flexmix_test <- computeSNR(processed, update_correspondence = TRUE, package = "flexmix")
  expect_equal(flexmix_test$CorrespondenceMatrix$snr, c(
    1, 1.46662215920315, 0, 0, 2.00872721014486, 0, 3.08641724395242,
    3.80602499772131, 3.31489921308867, 15.2477686105868, 2.40921027994661,
    1.41702657054023, 2.9132994473942
  ))
})
