test_that("generateMetapeaks works",{

  # get input data
  path <- system.file("extdata/Example_data.imzML", package = "gutenTAG")
  raw <- readMSIData(path)
  pre <- preProcess(raw, cores = 2)
  peaks <- peakDetection(pre, core = 2)

  # test that the object after peakDetection is still the correct object class
  expect_s4_class(peaks, "MSProcessedImagingExperiment")

  # generate metapeaks
  cur_test <- generateMetapeaks(peaks)

  # test that metapeak center is correct
  expect_equal(cur_test$metapeaks$center[1:10], c(903.138189115978, 906.171296812996, 911.900663726184, 915.997366639989,
                                                   919.190190221908, 921.226647010216, 923.182677484328, 925.219974994659,
                                                   927.994573754, 936.123016678659))

  # test that metapeak max is correct
  expect_equal(cur_test$metapeaks$max[1:10], c(903.241638183594, 905.2470703125, 912.266174316406, 914.77294921875,
                                               919.173828125, 921.234985351562, 923.184692382812, 925.245849609375,
                                               927.195617675781, 937.779907226562))

  # test that metapeak width is correct
  expect_equal(cur_test$metapeaks$width[1:10], c(0.167565628749472, 0.117917435010306, 0.263757941662039, 0.176874452584152,
                                               0.0372352953387212, 0.0713699159742911, 0.0434434444702502, 0.0093088238346803,
                                               0.117917435010306, 0.39408827507279))


  # test that max value of propagation selection is correct
  expect_equal(max(cur_test$propagation_selection), 158)

  # test if character in first argument throws error
  expect_error(generateMetapeaks("test"))

  # check if fixed limits works
  limits <- 0.5
  fixed_lim_test <- generateMetapeaks(x = peaks, fixed.limits = limits)

  # check the metapeak widths are all the same
  expect_true(all(fixed_lim_test$metapeaks$width == fixed_lim_test$metapeaks$width[1]))
  # check that the metapeak limits range is double the user-defined limits
  expect_equal(mean(fixed_lim_test$metapeaks$limits[, 2] - fixed_lim_test$metapeaks$limits[, 1]), limits*2)

})
