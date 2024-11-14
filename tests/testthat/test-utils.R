# .cleanPanel
test_that(".cleanPanel works", {
  panel_path <- system.file("extdata/ref_list.csv", package = "gutenTAG")
  panel <- read.delim(panel_path, sep = ",", header = TRUE, col.names = c("Name", "FeatureMass"))
  panel <- dplyr::arrange(panel, panel$FeatureMass)

  cur_test <- .cleanPanel(panel)


  # test that processed data is a list
  expect_type(cur_test, "list")

  # test that dimensions are corrent
  expect_equal(dim(cur_test), c(13, 2))

  # test that order of names is correct
  expect_equal(cur_test$Name, c(
    "CD98", "NFKB", "FN1", "CD73", "beta.actin", "VIM", "Collagen.1A1", "AASM", "Caveolin1",
    "NapsinA", "CK7", "ATP5a", "HLA.ABC"
  ))

  # test that image values are correct
  expect_equal(cur_test$FeatureMass, c(
    943.59, 974.53, 1068.6, 1088.57, 1206.72, 1230.84, 1234.87,
    1251.68, 1356.68, 1386.72, 1432.77, 1564.76, 1569.8
  ))

  # test if character in first argument throws error
  expect_error(.cleanPanel("panel_path"))
})

# .correctCoordinates. doesn't add coverage
test_that(".correctCoordinates works", {
  path <- system.file("extdata/Example_data.imzML", package = "gutenTAG")
  panel_path <- system.file("extdata/ref_list.csv", package = "gutenTAG")
  panel <- readPanel(path = panel_path)
  raw <- readMSIData(path)
  pre <- preProcess(raw, cores = 2)

  coords <- as.data.frame(Cardinal::coord(pre))

  cur_test <- .correctCoordinates(coords)

  # test that processed data is a list
  expect_type(cur_test, "list")

  # test that dimensions are corrent
  expect_equal(dim(cur_test), c(256L, 2L))

  # test that order of names is correct
  expect_equal(head(cur_test$x, 10), c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1))
  expect_equal(head(cur_test$y, 10), c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10))

  # test if character in first argument throws error
  expect_error(.cleanPanel("panel_path"))
})

#
