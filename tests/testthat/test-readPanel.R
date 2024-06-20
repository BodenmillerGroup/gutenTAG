test_that("readPanel works",{

  # get input data
  path <- system.file("extdata/Example_data.imzML", package = "gutenTAG")
  panel_path <- system.file("extdata/ref_list.csv", package = "gutenTAG")

  expect_silent(cur_test <- readPanel(path = panel_path))

  # test that output is correct is correct
  expect_equal(cur_test[[2]], c(943.59, 974.53, 1068.6, 1088.57, 1206.72, 1230.84, 1234.87,
                                1251.68, 1356.68, 1386.72, 1432.77, 1564.76, 1569.8))
  expect_equal(cur_test[[1]], c("CD98", "NFKB", "FN1", "CD73", "beta.actin", "VIM", "Collagen.1A1",
                                "AASM", "Caveolin1", "NapsinA", "CK7", "ATP5a", "HLA.ABC"))

  # test that dimensions are correct
  expect_equal(dim(cur_test), c(13, 2))

  # test if input is a string
  not_a_string <- 100L
  expect_error(readPanel(not_a_string))

})

