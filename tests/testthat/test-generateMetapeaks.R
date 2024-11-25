test_that("generateMetapeaks works",{

  # get input data
  path <- system.file("extdata/Example_data.imzML", package = "gutenTAG")
  raw <- readMSIData(path)
  pre <- preProcess(raw, cores = 2)
  peaks <- peakDetection(pre, core = 2)

  # test that the object after peakDetection is still the correct object class
  expect_s4_class(peaks, "MSImagingExperiment")

  # generate metapeaks
  cur_test <- generateMetapeaks(peaks, hist_smooth_factor = 1)

  # test that metapeak center is correct
  expect_equal(cur_test$metapeaks$center[1:10], c(903.09, 906.31, 913.88, 917.16,
                                                  920.11, 923.15, 926.09, 932.97,
                                                  937.78, 941.23), tolerance = 0.1)

  # test that metapeak max is correct
  expect_equal(cur_test$metapeaks$max[1:10], c(903.19, 905.19, 914.72, 917.17,
                                               919.17, 923.13, 925.19, 929.20,
                                               937.78, 941.18), tolerance = 0.1)

  # test that metapeak width is correct
  expect_equal(cur_test$metapeaks$width[1:10], c(0.338231932202095, 0.124122184279222, 0.235831470157998, 0.148947981080111,
                                                 0.117917435010306, 0.449944617943485, 0.232730795454847, 0.0124094985378317,
                                                 0.40339709890747, 0.505797560951566), tolerance = 0.01)


  # test that max value of propagation selection is correct
  expect_equal(max(cur_test$propagation_selection), 245)

  # test if character in first argument throws error
  expect_error(generateMetapeaks("test"))

  # check if fixed limits works
  limits <- 0.5
  fixed_lim_test <- generateMetapeaks(x = peaks, fixed.limits = limits, hist_smooth_factor = 1)

  # check the metapeak widths are all the same
  expect_true(all(fixed_lim_test$metapeaks$width == fixed_lim_test$metapeaks$width[1]))
  # check that the metapeak limits range is double the user-defined limits
  expect_equal(mean(fixed_lim_test$metapeaks$limits[, 2] - fixed_lim_test$metapeaks$limits[, 1]), limits*2)

  # check that it doesn't work when fixed.limits isn't a number
  expect_error(generateMetapeaks(x = peaks, fixed.limits = "limits"))

  # check that logical input produces error
  expect_error(generateMetapeaks(x = peaks, fixed.limits = TRUE))


})
