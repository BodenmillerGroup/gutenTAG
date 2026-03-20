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
  expect_equal(cur_test@listData$test@.Data[1:5,1:5,1], structure(c(17.4220487250521, 12.3030317972148, 11.022383809441, 
                                                                    14.4095816115594, 15.3221911771007, 15.5347563384415, 17.6477765907914, 
                                                                    15.8208378409457, 12.250847244694, 16.1245279689169, 16.6012484948729, 
                                                                    15.6629526063614, 20.338495529729, 18.2699194271189, 13.4859455343943, 
                                                                    13.599311601724, 12.8299417579049, 13.597775511662, 20.0133271685664, 
                                                                    14.8645459066664, 15.4860512443758, 18.3871916891564, 14.9019565797177, 
                                                                    15.8036435582998, 16.776347033104), dim = c(5L, 5L)),
               tolerance = 0.001)

  # test if character in first argument throws error
  expect_error(asCytoImageList("test"))


})
