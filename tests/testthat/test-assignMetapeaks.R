test_that("assignMetapeaks works",{

  # get input data
  path <- system.file("extdata/Example_data.imzML", package = "gutenTAG")
  panel_path <- system.file("extdata/ref_list.csv", package = "gutenTAG")
  panel <- readPanel(path = panel_path)
  raw <- readMSIData(path)
  pre <- preProcess(raw, cores = 2)
  peaks <- peakDetection(pre, core = 2)
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
  expect_equal(dim(cur_test$Untargeted$UntargetedCorrespondence), c(158, 7))
  expect_equal(dim(cur_test$Untargeted$UntargetedIntensity), c(256, 158))
  # check that number of untargeted metapeaks is same as in watershed output
  expect_equal(dim(cur_test$Untargeted$UntargetedIntensity)[2], max(metapeaks$propagation_selection))

  # check if output of intensity df is the same
  expect_equal(cur_test$IntensityDF[1][1:20,], c(7.82970978419113, 7.46634943048983, 7.32187611962462, 6.54747422274114,
                                                7.49722025536448, 3.55251445305979, 5.75086256527666, 3.51621186218398,
                                                5.27250246930296, 3.86230533109536, 3.44947707629045, 9.10591335737913,
                                                8.21399509393651, 7.16647250531501, 6.75161308428109, 7.91584977381878,
                                                6.24461353570202, 9.09662347324635, 6.15357788356246, 6.6490763080287))

  ### Uncomment when CorrespondenceMatrix$expected_mz_location is correctly changed from character to numeric ###

  # check if output of correspondence matrix is the same
  #expect_equal(cur_test$CorrespondenceMatrix$mean, c(7.84679001732988, 16.1002595411377, NA, NA, 26.0409182022813,
  #                                                   NA, 14.4482423456346, 9.41933313316087, 4.88193039649564, 14.7482120965546,
  #                                                   1.71918838044124, 6.78957980506437, 8.7250571739829))

  # check if spatial coordinates are correct
  expect_equal(range(cur_test$SpatialCoords), c(1, 16))


  # test if character in arguments throws error
  expect_error(assignMetapeaks("test", pre = pre, refList = panel))
  expect_error(assignMetapeaks(metapeaks, pre = "test", refList = panel))
  expect_error(assignMetapeaks(metapeaks, pre = pre, refList = "panel"))

  # test that mz_threshold should be a single numeric
  expect_error(assignMetapeaks(metapeaks, pre = pre, refList = panel, mz_threshold = "test"))


})
