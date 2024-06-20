test_that("asMSImagingExperiment works",{

  # get input data
  path <- system.file("extdata/Example_data.imzML", package = "gutenTAG")
  panel_path <- system.file("extdata/ref_list.csv", package = "gutenTAG")
  panel <- readPanel(path = panel_path)
  raw <- readMSIData(path)
  pre <- preProcess(raw, cores = 2)
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
  expect_equal(round(fData(cur_test)[[2]], 2), c(943.59, 974.53, 1068.60, 1088.57, 1206.72, 1230.84, 1234.87, 1251.68, 1356.68, 1386.72,
                                      1432.77, 1564.76, 1569.8))

  # test that values in intensity are correct
  expect_equal(iData(cur_test)[1:5, 1:5], structure(c(7.82970978419113, 17.7493402240123, 38.21930389974,
                                                      38.121208173945, 14.2061520368805, 7.46634943048983, 15.5531491674357,
                                                      41.2465855848849, 31.0258658353397, 10.2214001836892, 7.32187611962462,
                                                      15.5819123732679, 44.6965441458399, 32.5160404293837, 11.5167339794897,
                                                      6.54747422274114, 16.5104467390468, 28.7510100965559, 15.9670806372939,
                                                      4.01539500585915, 7.49722025536448, 18.7082550227178, 25.2541078944199,
                                                      7.45608285560947, 6.49318138965854), dim = c(5L, 5L), dimnames = list(
                                                        c("CD98", "NFKB", "beta.actin", "Collagen.1A1", "AASM"),
                                                        NULL)))



  # test if character in first argument throws error
  expect_error(asMSImagingExperiment("test"))


})
