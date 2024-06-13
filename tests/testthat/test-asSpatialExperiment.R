test_that("asSpatialExperiment works",{

  # get input data
  path <- system.file("extdata/Example_data.imzML", package = "gutenTAG")
  panel_path <- system.file("extdata/ref_list.csv", package = "gutenTAG")
  panel <- readPanel(path = panel_path)
  raw <- readMSIData(path)
  pre <- preProcess(raw, cores = 2)
  peaks <- peakDetection(pre, core = 2)
  metapeaks <- generateMetapeaks(peaks)
  processed <- assignMetapeaks(metapeaks, pre = pre, refList = panel)

  # test that processed data is a list
  expect_type(processed, "list")

  # make SpatialExperiment object
  cur_test <- asSpatialExperiment(processed)

  # test that type is SpatialExperiment
  expect_s4_class(cur_test, "SpatialExperiment")


  ### Uncomment when CorrespondenceMatrix$expected_mz_location is correctly changed from character to numeric ###

  # test that intensity df is the same
  #expect_equal(cur_test@assays@data$intenity[1:20], c(7.82970978419113, 17.7493402240123, 0, 0, 38.21930389974, 0,
  #                                                    38.121208173945, 14.2061520368805, 7.78127437244719, 10.8479657773842,
  #                                                    1.55775389316223, 9.09659022610243, 14.9059079381073, 7.46634943048983,
  #                                                    15.5531491674357, 0, 0, 41.2465855848849, 0, 31.0258658353397))

  # test that first entries into cur_test are equal to first 20 entries into processed intensity df
  expect_equal(cur_test@assays@data$intensity[1, ][1:20], processed$IntensityDF[[1]][1:20])

  # test if character in first argument throws error
  expect_error(asSpatialExperiment("test"))


})
