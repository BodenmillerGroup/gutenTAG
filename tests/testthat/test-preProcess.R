test_that("preProcess works",{

  path <- system.file("extdata/Example_data/Example_data.imzML", package = "gutenTAG")
  raw <- readMSIData(path)

  # test that the object loaded in is of the correct class
  expect_s4_class(raw, "MSImagingExperiment")

  # test that message is generated during preProcessing
  #expect_message(cur_test <- preProcess(raw, cores = 2))
  cur_test <- preProcess(raw)

  # test that the object after preProcess is still the correct object class
  expect_s4_class(cur_test, "MSImagingExperiment")

  # test first 20 values of intensity df
  expect_equal(Cardinal::spectra(cur_test)[1:20],
               c(0, 0.0344849614269698, 0.017047294500572, 0, 0, 0.0239869996341981,
                 0.0554512608941657, 0.119392982958655, 0.190402790127705, 0.161389544976431,
                 0.068514460283599, 0, 0.130981346146514, 0.433087156557264, 0.851486432625933,
                 1.15127021255953, 1.15668061083919, 0.962086743849345, 0.6767970548317,
                 0.329735412157482),
               tolerance = 0.001)

  # test if character in first argument throws error
  expect_error(preProcess("test"))

  # test if invalid BPPARAM throws error
  expect_error(preProcess(raw, BPPARAM = "not_a_param"))

})

