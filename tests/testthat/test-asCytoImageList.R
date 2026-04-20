test_that("asCytoImageList works",{

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
  expect_equal(cur_test@listData$test@.Data[1:5,1:5,1], structure(c(18.60279186794, 13.3238938839062, 12.1614539299946, 
                                                                    14.954563035264, 16.4515142167085, 16.0857977316972, 18.007727014023, 
                                                                    16.6355121776616, 12.9747603696865, 16.7355243061082, 17.9647605811798, 
                                                                    16.5496371996611, 22.6771209023691, 19.8142685150797, 14.9410920042159, 
                                                                    15.0788161493883, 13.3184652430064, 14.8716370690039, 21.8627752504725, 
                                                                    15.2202371885168, 16.2613473790588, 19.0847197894317, 16.0117522631747, 
                                                                    16.6318727212659, 17.627609189568), dim = c(5L, 5L)),
               tolerance = 0.001)

  # test if character in first argument throws error
  expect_error(asCytoImageList("test"))


})
