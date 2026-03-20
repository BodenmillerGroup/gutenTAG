test_that("generateMetapeaks works",{

  # get input data
  path <- system.file("extdata/Example_data.imzML", package = "gutenTAG")
  raw <- readMSIData(path)
  pre <- preProcess(raw)
  peaks <- peakDetection(pre)

  # test that the object after peakDetection is still the correct object class
  expect_s4_class(peaks, "MSImagingExperiment")

  # generate metapeaks
  cur_test <- generateMetapeaks(peaks, hist_smooth_factor = 1, sparsity = 3)

  # test that metapeak center is correct
  expect_equal(cur_test$metapeaks$center[1:10], c(903.09, 906.31, 913.88, 917.16,
                                                  920.11, 923.15, 926.09, 932.97,
                                                  937.78, 941.23), tolerance = 0.1)

  # test that metapeak max is correct
  expect_equal(cur_test$metapeaks$max[1:10], c(903.19, 905.19, 914.72, 917.17,
                                               919.17, 923.13, 925.19, 929.20,
                                               937.78, 941.18), tolerance = 0.1)

  # test that metapeak width is correct
  expect_equal(cur_test$metapeaks$width[1:10], c(0.338231932202095, 0.366158403706136, 0.0557033490531052, 0.148944581217498, 
                                                 0.0527522683049305, 0.00620474926891583, 0.117917435010306, 0.32581903380165, 
                                                 0.0124128984004448, 0.0155135731035961), tolerance = 0.01)


  # test that max value of propagation selection is correct
  expect_equal(max(cur_test$propagation_selection), 225)

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

  # check that sparsity param fails if not a positive number
  expect_error(generateMetapeaks(x = peaks, fixed.limits = TRUE, sparsity = TRUE))
  expect_error(generateMetapeaks(x = peaks, fixed.limits = TRUE, sparsity = "2"))
  expect_error(generateMetapeaks(x = peaks, fixed.limits = TRUE, sparsity = as.factor(2)))
  expect_error(generateMetapeaks(x = peaks, fixed.limits = TRUE, sparsity = -1))

  # check sparsity output is expected
  sparsity_test <- generateMetapeaks(x = peaks, fixed.limits = limits, hist_smooth_factor = 1, sparsity = 1)
  expect_equal(sparsity_test$metapeaks$max[1:10] , c(903.185974121094, 905.19140625, 907.196838378906, 912.266174316406,
                   914.71728515625, 917.168395996094, 919.173828125, 921.234985351562,
                   923.129028320312, 925.190185546875), tolerance = 0.1)
  # (note: sparsity = 1 → lower value → more closely-spaced metapeaks)




})
