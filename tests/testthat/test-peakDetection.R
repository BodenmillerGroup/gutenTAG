test_that("peakDetection works",{

  # get input data
  path <- system.file("extdata/Example_data.imzML", package = "gutenTAG")
  raw <- Cardinal::readMSIData(path)
  pre <- preProcess(raw)

  # test that the object loaded in is of the correct class
  expect_s4_class(pre, "MSImagingExperiment")
  expect_false(Cardinal::isCentroided(pre))

  # test that message is generated during peakDetection
  cur_test <- peakDetection(raw, snr = 3)

  # test that the object after peakDetection is still the correct object class
  expect_s4_class(cur_test, "MSImagingExperiment")
  expect_false(Cardinal::isCentroided(pre))

  # test that dimensions of intensity after peakDetection are the same
  expect_equal(dim(ProtGenerics::spectra(cur_test)), c(15800, 64))

  expect_equal(
    ProtGenerics::spectra(cur_test)[5510, 1:10],
    c(9.19874954223633, 9.0649995803833, 0, 0, 5.33249998092651, 
      2.2574999332428, 2.23874998092651, 3.65000009536743, 5.74625015258789, 
      8.36999988555908),
    tolerance = 0.001
  )

  # test that coords are correct
  expect_equal(dim(Cardinal::coord(cur_test)), c(64, 2))

  # test if errors in arguments throws errors
  expect_error(peakDetection(x = "test", snr = 3, win = 50))
  expect_error(peakDetection(x = test, snr = -1, win = 50))
  expect_error(peakDetection(x = cur_test, snr = 3, win = 50, BPPARAM = "not_a_param"))


})
