test_that("peakDetection works",{

  # get input data
  path <- system.file("extdata/Example_data.imzML", package = "maldipackage")
  raw <- readMSIData(path)
  pre <- preProcess(raw, cores = 2)

  # test that the object loaded in is of the correct class
  expect_s4_class(pre, "MSContinuousImagingExperiment")

  # test that message is generated during peakDetection
  #expect_message(cur_test <- peakDetection(raw, snr = 3, cores = 2))
  cur_test <- peakDetection(raw, snr = 3, cores = 2)

  # test that the object after peakDetection is still the correct object class
  expect_s4_class(cur_test, "MSProcessedImagingExperiment")

  # test that pre-processing operations are correct
  expect_equal(names(cur_test@processing), "peakPick")

  # test that dimensions of intensity after peakDetection are the same
  expect_equal(dim(cur_test@imageData$data$intensity), c(15800, 256))

  # test if character in first argument throws error
  expect_error(peakDetection("test", cores = 2))


})
