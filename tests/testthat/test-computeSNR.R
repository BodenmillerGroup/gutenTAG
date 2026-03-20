test_that("computeSNR works",{

  # get input data
  path <- system.file("extdata/Example_data.imzML", package = "gutenTAG")
  panel_path <- system.file("extdata/ref_list.csv", package = "gutenTAG")
  panel <- readPanel(path = panel_path)
  raw <- readMSIData(path)
  pre <- preProcess(raw)
  peaks <- peakDetection(pre)
  metapeaks <- generateMetapeaks(peaks)
  processed <- assignMetapeaks(metapeaks, pre = pre, refList = panel)

  # compute Geary's C
  cur_test <- computeSNR(processed, update_correspondence = TRUE, package = "mclust")

  # test if character in first argument throws error
  expect_error(computeSNR("processed"))

  # Test correspondence matrix
  expect_true("snr" %in% colnames(cur_test$CorrespondenceMatrix))
  expect_equal(cur_test$CorrespondenceMatrix$snr, c(1.46674417788158, 1.68252354822219, 1.89190839866574, 0, 3.99920021456553, 
                                                    0, 2.44622835834056, 2.90072853664465, 3.81615438430949, 0, 1.54626104953847, 
                                                    1.38225470794779, 3.10914667005805))

  # test flexmix snr
  set.seed(123)
  flexmix_test <- computeSNR(processed, update_correspondence = TRUE, package = "flexmix")
  expect_equal(flexmix_test$CorrespondenceMatrix$snr, c(1, 1.65067019405939, 1, 0, 1, 0, 2.51257502580874, 1, 4.25386455758169, 
                                                        0, 1.29670203517052, 1.36948499744764, 3.10914667005805), tolerance = 0.1)

})

