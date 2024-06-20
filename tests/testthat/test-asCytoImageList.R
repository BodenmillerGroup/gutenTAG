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
  expect_equal(cur_test@listData$test@.Data[1:5,1:5,1], structure(c(7.82970978419113, 6.24461353570202, 4.99465676520053,
                                                                    7.6344321582305, 8.274748861544, 7.46634943048983, 9.09662347324635,
                                                                    6.80983928972767, 5.37533742996382, 8.04719719339299, 7.32187611962462,
                                                                    6.15357788356246, 6.29577848525015, 8.4844615690162, 5.64467915346895,
                                                                    6.54747422274114, 6.6490763080287, 6.72442700638513, 11.8954635507536,
                                                                    9.27353112319082, 7.49722025536448, 9.19165723131113, 7.15812443090326,
                                                                    7.50190937599801, 7.66192188916043), dim = c(5L, 5L)))

  # test if character in first argument throws error
  expect_error(asCytoImageList("test"))


})
