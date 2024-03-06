test_that("pca works",{

  # get input data
  path <- system.file("extdata/Example_data.imzML", package = "maldipackage")
  panel_path <- system.file("extdata/ref_list.csv", package = "maldipackage")
  panel <- readPanel(path = panel_path)
  raw <- readMSIData(path)
  pre <- preProcess(raw, cores = 2)
  peaks <- peakDetection(pre, core = 2)
  metapeaks <- generateMetapeaks(peaks)
  processed <- assignMetapeaks(metapeaks, pre = pre, refList = panel)

  # test that processed data is a list
  expect_type(processed, "list")

  # test that entries into intensity dataframe are correct
  expect_equal(processed$IntensityDF[1:20, ][[1]], c(7.82970978419113, 7.46634943048983, 7.32187611962462, 6.54747422274114,
                                                     7.49722025536448, 3.55251445305979, 5.75086256527666, 3.51621186218398,
                                                     5.27250246930296, 3.86230533109536, 3.44947707629045, 9.10591335737913,
                                                     8.21399509393651, 7.16647250531501, 6.75161308428109, 7.91584977381878,
                                                     6.24461353570202, 9.09662347324635, 6.15357788356246, 6.6490763080287))

  ### Uncomment when CorrespondenceMatrix$expected_mz_location is correctly changed from character to numeric ###

  # make CytoImageList object
  #cur_test <- computePCA(processed$IntensityDF[-c(3,4,6)], comp = 2, seed = 123, scree = FALSE)

  # test that class of s3 object is correct
  #expect_equal(class(cur_test), c("irlba_prcomp", "prcomp"))
  #expect_s3_class(cur_test, c("irlba_prcomp", "prcomp"))
#
  ## test that loadings are the same
  #expect_equal(cur_test$rotation[, 1], c(0.0398562503389398, 0.338502919521346, 0.416082164272529, 0.396842233137603,
  #                                       0.40498131672024, 0.417928083264583, -0.0572137926221178, 0.0227850771200478,
  #                                       0.233816269605792, 0.395073096964957))
#
  ## test that sdev is the same
  #expect_equal(cur_test$sdev, c(2.19855868374997, 1.30258144129773, 1.02141994170742))
#
#
  ## test if arguments are in correct format
  #expect_error(computePCA("test", comp = 5, seed = 123, scree = FALSE))
  ## test if breaks when just input processed
  #expect_error(computePCA(processed, comp = 5, seed = 123, scree = FALSE))
  #expect_error(computePCA(processed$IntensityDF, comp = "5", seed = 123, scree = FALSE))

  # expect warning in normal case
  #expect_warning(computePCA(processed$IntensityDF, comp = 3, seed = 123, scree = FALSE))

})
