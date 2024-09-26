test_that("asCytoImageList works",{

  # am I allowed to do this? ask Nils
  library(cytomapper)

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

  # make CytoImageList object
  cur_test <- asCytoImageList(processed)

  # test that type is CytoImageList
  expect_s4_class(cur_test, "CytoImageList")

  # test that channelNames are correct
  expect_equal(length(channelNames(cur_test)), 13)

  ### Uncomment when introduced change that columns are ordered correctly (in m/z order as in panel) ###

  #expect_equal(channelNames(cur_test), c("CD98", "NFKB", "FN1", "CD73", "beta.actin", "VIM", "Collagen.1A1", "AASM", "Caveolin1",
  #                                       "NapsinA", "CK7", "ATP5a", "HLA.ABC"))

  # test that image values are correct
  expect_equal(cur_test@listData$test@.Data[1:5,1:5,1], structure(c(9.78751411382999, 7.49889435489177, 5.89596411683664,
                                                                    9.40937495640716, 9.43915222733905, 8.83914708481089, 11.2326644673979,
                                                                    8.31376932515358, 6.14837567470535, 9.71965830322728, 9.20862181033298,
                                                                    7.62966172969292, 9.21913374175065, 10.6457146936996, 6.56208110350116,
                                                                    7.82156800846123, 7.37275919838429, 8.28032021970029, 13.6156501907851,
                                                                    8.99621765710804, 8.51768733290817, 10.8724118537528, 8.01363079998227,
                                                                    8.40351981435345, 9.27146306773115), dim = c(5L, 5L)),
               tolerance = 0.001)

  # test if character in first argument throws error
  expect_error(asCytoImageList("test"))


})
