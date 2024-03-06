test_that("scaleData works",{

  # get input data
  path <- system.file("extdata/Example_data.imzML", package = "maldipackage")
  panel_path <- system.file("extdata/ref_list.csv", package = "maldipackage")
  panel <- readPanel(path = panel_path)
  raw <- readMSIData(path)
  pre <- preProcess(raw, cores = 2)
  peaks <- peakDetection(pre, core = 2)
  metapeaks <- generateMetapeaks(peaks)
  processed <- assignMetapeaks(metapeaks, pre = pre, refList = panel)


  ### Uncomment when CorrespondenceMatrix$expected_mz_location is correctly changed from character to numeric ###


  # test if corrected standard deviation scaling works
  #cur_test <- scaleData(x = processed, method = "corsd")
#
  ## test that output is correct
  #expect_equal(cur_test$Collagen.1A1[1:10], c(5.15714027978038, 4.19726315295148, 4.3988580076524, 2.16006990988733,
  #                                            1.00867908090921, 1.06323103020873, 1.08744979291011, 1.09104548141196,
  #                                            1.09582715165859, 0.790123644452622))
#
  ## test that dimensions of output is correct
  #expect_equal(dim(cur_test), c(256L, 13L))
#
#
  ## test if geary scaling works
  #cur_test <- scaleData(x = processed, method = "geary")
#
  ## test that output is correct
  #expect_equal(cur_test$Collagen.1A1[1:10], c(2.85579363046194, 2.32425660875796, 2.43589082282927, 1.19615010555257,
  #                                            0.558561361174259, 0.588769790824634, 0.602181058408614, 0.604172189881416,
  #                                            0.606820065000659, 0.437535135499798))
#
  ## test that dimensions of output is correct
  #expect_equal(dim(cur_test), c(256L, 13L))
#
#
  ## test if character in first argument throws error
  #expect_error(scaleData(x = "processed", method = "corsd"))


})
