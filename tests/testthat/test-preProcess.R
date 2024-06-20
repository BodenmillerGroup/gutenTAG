test_that("preProcess works",{

  path <- system.file("extdata/Example_data.imzML", package = "gutenTAG")
  raw <- readMSIData(path)

  # test that the object loaded in is of the correct class
  expect_s4_class(raw, "MSContinuousImagingExperiment")

  # test that message is generated during preProcessing
  #expect_message(cur_test <- preProcess(raw, cores = 2))
  cur_test <- preProcess(raw, cores = 2)

  # test that the object after preProcess is still the correct object class
  expect_s4_class(cur_test, "MSContinuousImagingExperiment")

  # test that pre-processing operations are correct
  expect_equal(names(cur_test@processing), c("normalize", "smoothSignal", "reduceBaseline"))

  # test first 20 values of intensity df
  expect_equal(cur_test@imageData$data$intensity[1:20], c(0.188345944717598, 0.161569424639892, 0.0967960608948741, 0.0245269394553258,
                                                          0, 0.0018638128854811, 0.0473214395395325, 0.0986611763036347,
                                                          0.14955623441025, 0.115596724204139, 0.0435247842738429, 0, 0.125678415781081,
                                                          0.418871193383176, 0.799962683685664, 1.05650252611367, 1.09211283819106,
                                                          0.915272796475069, 0.62950584957254, 0.284536804237329))

  # test if character in first argument throws error
  expect_error(preProcess("test", cores = 2))

})

