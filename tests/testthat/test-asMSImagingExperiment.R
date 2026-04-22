test_that("asMSImagingExperiment works",{

  # get input data
  path <- system.file("extdata/Example_data.imzML", package = "gutenTAG")
  panel_path <- system.file("extdata/ref_list.csv", package = "gutenTAG")
  panel <- readPanel(path = panel_path)
  raw <- readMSIData(path)
  pre <- preProcess(raw)
  peaks <- peakDetection(pre)
  metapeaks <- generateMetapeaks(peaks)
  processed <- assignMetapeaks(metapeaks, pre = pre, refList = panel)

  # test that processed data is a list
  expect_type(processed, "list")

  # make MSImagingExperiment object
  cur_test <- asMSImagingExperiment(processed)


  # test that type is MSImagingExperiment
  expect_s4_class(cur_test, "MSImagingExperiment")

  # test that mz values are the same
  expect_equal(ProtGenerics::mz(cur_test),
               c(943.59, 974.53, 1068.6, 1088.57, 1206.72, 1230.84, 1234.87,
                 1251.68, 1356.68, 1386.72, 1432.77, 1564.76, 1569.8),
               tolerance = 0.01)

  # test that values in intensity are correct
  expect_equal(ProtGenerics::spectra(cur_test)[1:5, 1:5],
               structure(c(18.60279186794, 11.864188202163, 8.84359502901302, 
                           0, 124.052120627902, 16.0857977316972, 10.8452486371537, 8.49246563090683, 
                           0, 124.132562738403, 17.9647605811798, 9.34187594768916, 9.15812074103582, 
                           0, 138.737857074253, 15.0788161493883, 12.7694536418265, 7.52165699686672, 
                           0, 99.2124832195187, 16.2613473790588, 11.1005731458917, 10.2079352756953, 
                           0, 70.659911492863), dim = c(5L, 5L), dimnames = list(c("CD98", 
                                                                                   "NFKB", "FN1", "CD73", "beta.actin"), NULL)),
               tolerance = 0.01)

  ## test remove NA clause
  new_test <- asMSImagingExperiment(processed, remove.na = TRUE)

  # check shape is correct
  expect_equal(as.numeric(dim(new_test)), c(10, 64))

  # check there are no more zero columns
  expect_true(all(rowSums(ProtGenerics::spectra(new_test)) > 0))

  # check NA values were correctly removed
  tmp <- fData(cur_test)
  expect_equal(fData(new_test), fData(cur_test)[!is.na(tmp$observed_mz), ])

  # test if character in first argument throws error
  expect_error(asMSImagingExperiment("test"))


})
