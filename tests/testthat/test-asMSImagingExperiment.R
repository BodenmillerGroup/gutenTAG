test_that("asMSImagingExperiment works",{

  # get input data
  path <- system.file("extdata/Example_data.imzML", package = "gutenTAG")
  panel_path <- system.file("extdata/ref_list.csv", package = "gutenTAG")
  panel <- readPanel(path = panel_path)
  raw <- readMSIData(path)
  pre <- preProcess(raw)
  peaks <- peakDetection(pre, core = 2)
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
               structure(c(9.78751411382999, 18.9425446633928, 0, 0, 50.5297566319526,
                           8.83914708481089, 16.8270374577648, 0, 0, 52.6146794330136, 9.20862181033298,
                           16.8391295927981, 0, 0, 59.2717356611516, 7.82156800846123, 19.7860939843446,
                           0, 0, 38.417457849565, 8.51768733290817, 20.4836879757656, 0,
                           0, 32.0163733894236), dim = c(5L, 5L), dimnames = list(c("CD98",
                                                                                    "NFKB", "FN1", "CD73", "beta.actin"), NULL)),
               tolerance = 0.01)

  ## test remove NA clause
  new_test <- asMSImagingExperiment(processed, remove.na = TRUE)

  # check shape is correct
  expect_equal(as.numeric(dim(new_test)), c(10, 256))

  # check there are no more zero columns
  expect_true(all(rowSums(ProtGenerics::spectra(new_test)) > 0))

  # check NA values were correctly removed
  tmp <- fData(cur_test)
  expect_equal(fData(new_test), fData(cur_test)[!is.na(tmp$observed_mz), ])

  # test if character in first argument throws error
  expect_error(asMSImagingExperiment("test"))


})
