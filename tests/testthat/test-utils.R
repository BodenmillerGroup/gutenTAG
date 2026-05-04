# .cleanPanel
test_that(".cleanPanel works",{

  panel_path <- system.file("extdata/ref_list.csv", package = "gutenTAG")
  panel <- read.delim(panel_path, sep = ",", header = TRUE, col.names = c("Name", "FeatureMass"))
  panel <- dplyr::arrange(panel, panel$FeatureMass)

  cur_test <- .cleanPanel(panel)


  # test that processed data is a list
  expect_type(cur_test, "list")

  # test that dimensions are corrent
  expect_equal(dim(cur_test), c(13, 3))
  expect_true("OriginalName" %in% colnames(cur_test))

  # test that order of names is correct
  expect_equal(cur_test$Name, c("CD98", "NFKB", "FN1", "CD73", "beta.actin", "VIM", "Collagen.1A1", "AASM", "Caveolin1",
                           "NapsinA", "CK7", "ATP5a", "HLA.ABC"))

  # test that image values are correct
  expect_equal(cur_test$FeatureMass, c(943.59, 974.53, 1068.6, 1088.57, 1206.72, 1230.84, 1234.87,
                                       1251.68, 1356.68, 1386.72, 1432.77, 1564.76, 1569.8))

  # test if character in first argument throws error
  expect_error(.cleanPanel("panel_path"))


})

test_that(".cleanPanel captures OriginalName before cleaning", {
  dirty_panel <- data.frame(
    Name        = c("beta-actin", "HLA/ABC", "CD 98"),
    FeatureMass = c(1206.72, 1569.8, 943.59)
  )
  result <- .cleanPanel(dirty_panel)
  # sorted by FeatureMass: CD 98 (943), beta-actin (1206), HLA/ABC (1569)
  expect_equal(result$OriginalName, c("CD 98", "beta-actin", "HLA/ABC"))
  expect_equal(result$Name,         c("CD98",  "beta.actin", "HLAABC"))
})

# .translate_coordinates doesn't add coverage
test_that(".translate_coordinates works",{

  path <- system.file("extdata/Example_data/Example_data.imzML", package = "gutenTAG")
  panel_path <- system.file("extdata/ref_list.csv", package = "gutenTAG")
  panel <- readPanel(path = panel_path)
  raw <- readMSIData(path)
  pre <- preProcess(raw)

  coords <- as.data.frame(Cardinal::coord(pre))

  cur_test <- .translate_coordinates(coords)

  # test that processed data is a list
  expect_type(cur_test, "list")

  # test that dimensions are corrent
  expect_equal(dim(cur_test), c(64L, 2L))

  # test that order of names is correct
  expect_equal(head(cur_test$x, 10), c(1, 1, 1, 1, 1, 1, 1, 1, 2, 2))
  expect_equal(head(cur_test$y, 10), c(1, 2, 3, 4, 5, 6, 7, 8, 1, 2))

  # test if character in first argument throws error
  expect_error(.cleanPanel("panel_path"))


})

# .choose_closest
test_that(".choose_closest works", {

  # picks the nearer candidate
  expect_equal(.choose_closest(c(98, 103), 100), 98)

  # candidates only to the right
  expect_equal(.choose_closest(c(102, 105), 100), 102)

  # candidates only to the left
  expect_equal(.choose_closest(c(95, 97), 100), 97)

  # tie — which.min returns first element
  expect_equal(.choose_closest(c(98, 102), 100), 98)

  # single candidate
  expect_equal(.choose_closest(105, 100), 105)

})

# .removeDuplicates
test_that(".removeDuplicates works", {

  # minimal fixture — two candidates for one expected mz (100) plus one
  # non-duplicate marker. The duplicate column uses the ".duplicate" suffix
  # matching real assignMetapeaks output. A non-duplicate column is required
  # because the function uses data.frame subsetting that drops dimensions when
  # only one column is kept.
  make_x <- function(mz_locations) {
    dup_cols  <- c(paste0("m", mz_locations[1]),
                   paste0("m", mz_locations[1], ".duplicate"))
    all_cols  <- c(dup_cols, "m200")
    list(
      IntensityDF = as.data.frame(
        setNames(replicate(3L, rnorm(10), simplify = FALSE), all_cols)
      ),
      CorrespondenceMatrix = data.frame(
        marker                = all_cols,
        expected_mz_location  = c(100, 100, 200),
        mz_location           = c(mz_locations, 200),
        stringsAsFactors      = FALSE
      ),
      FilteredDF = as.data.frame(
        setNames(replicate(3L, rnorm(10), simplify = FALSE), all_cols)
      )
    )
  }

  x_lr <- make_x(c(98, 103))   # one left, one right of expected mz 100

  # "closest" picks the nearer candidate (98, distance 2) over the farther (103, distance 3)
  result_closest <- .removeDuplicates(x_lr, shift = "closest")
  expect_equal(ncol(result_closest$IntensityDF), 2L)
  expect_equal(result_closest$CorrespondenceMatrix$mz_location, c(98, 200))

  # "closest" runs without a warning even when all candidates are on one side
  x_right_only <- make_x(c(102, 105))
  expect_silent(.removeDuplicates(x_right_only, shift = "closest"))

  x_left_only <- make_x(c(95, 97))
  expect_silent(.removeDuplicates(x_left_only, shift = "closest"))

  # "right" still selects the nearest right candidate
  result_right <- .removeDuplicates(x_lr, shift = "right")
  expect_equal(result_right$CorrespondenceMatrix$mz_location, c(103, 200))

  # "left" still selects the nearest left candidate
  result_left <- .removeDuplicates(x_lr, shift = "left")
  expect_equal(result_left$CorrespondenceMatrix$mz_location, c(98, 200))

  # "right" warns when no candidates exist to the right
  x_left_only2 <- make_x(c(95, 97))
  expect_warning(.removeDuplicates(x_left_only2, shift = "right"))

  # invalid shift errors
  expect_error(.removeDuplicates(x_lr, shift = "diagonal"))

})
