test_that("peakDetection works",{

  # get input data
  path <- system.file("extdata/Example_data.imzML", package = "maldipackage")
  raw <- readMSIData(path)
  pre <- preProcess(raw, cores = 2)

  # test that the object loaded in is of the correct class
  expect_s4_class(pre, "MSContinuousImagingExperiment")

  # test that message is generated during peakDetection
  cur_test <- peakDetection(raw, snr = 3, cores = 2)

  # test that the object after peakDetection is still the correct object class
  expect_s4_class(cur_test, "MSProcessedImagingExperiment")

  # test that pre-processing operations are correct
  expect_equal(names(cur_test@processing), "peakPick")

  # test that dimensions of intensity after peakDetection are the same
  expect_equal(dim(cur_test@imageData$data$intensity), c(15800, 256))

  expect_equal(cur_test@imageData$data$intensity@data[[1]], c(3.25, 2.45624995231628, 2.11999988555908, 1.58375000953674,
                                                              2.11249995231628, 2.29874992370605, 2.28125, 1.83124995231628,
                                                              2.37124991416931, 9.19874954223633, 2.18249988555908, 3.97874999046326,
                                                              1.64250004291534, 1.45000004768372, 3.42750000953674, 2.24125003814697,
                                                              1.45375001430511))

  # test that coords are correct
  expect_equal(dim(cur_test@elementMetadata@coord), c(256, 3))

  # test if errors in arguments throws errors
  expect_error(peakDetection(x = "test", snr = 3, win = 50, cores = 2))
  expect_error(peakDetection(x = test, snr = -1, win = 50, cores = 2))
  expect_error(peakDetection(x = test, snr = 3, win = 50, cores = 2))
  expect_error(peakDetection(x = test, snr = 3, win = 50, cores = -2))


})
