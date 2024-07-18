test_that("computeGearysC works",{

  # get input data
  path <- system.file("extdata/Example_data.imzML", package = "gutenTAG")
  panel_path <- system.file("extdata/ref_list.csv", package = "gutenTAG")
  panel <- readPanel(path = panel_path)
  raw <- readMSIData(path)
  pre <- preProcess(raw, cores = 2)
  peaks <- peakDetection(pre, core = 2)
  metapeaks <- generateMetapeaks(peaks)
  processed <- assignMetapeaks(metapeaks, pre = pre, refList = panel)


  # compute Geary's C
  expect_silent(cur_test <- computeGearysC(processed, verbose = FALSE))


  ### Uncomment when CorrespondenceMatrix$expected_mz_location is correctly changed from character to numeric ###

  # test that output is correct is correct
  #expect_equal(cur_test, c(0.710039425480746, 0.483790051545452, NaN, NaN, 0.435088747635843,
  #                         NaN, 0.434696662055006, 0.360060956047036, 0.301997161578108,
  #                         0.556339728384072, 0.759333786571679, 0.877932000280288, 0.391559739360737)
  #             )


  # test if character in first argument throws error
  expect_error(computeGearysC("processed", verbose = FALSE, update_correspondence = FALSE))
  expect_error(computeGearysC(processed, verbose = "FALSE", update_correspondence = FALSE))


  # Additional tests for higher coverage
  # Test with update_correspondence set to TRUE
  updated_processed <- computeGearysC(processed, verbose = FALSE, update_correspondence = TRUE)
  expect_true("GearysC" %in% colnames(updated_processed$CorrespondenceMatrix))

  # Test with missing values in the input data
  processed_with_na <- processed
  processed_with_na$IntensityDF[1, 1] <- NA
  expect_silent(cur_test_na <- computeGearysC(processed_with_na, verbose = FALSE))

  # Test with different k values in Knn function (modify the computeGearysC function temporarily for this test)
  #expect_silent({
  #  original_knn <- Knn
  #  assign("Knn", function(coords, k, verbose, indexType) original_knn(coords, k = 2, verbose = verbose, indexType = indexType), envir = .GlobalEnv)
  #  cur_test_k2 <- computeGearysC(processed, verbose = FALSE)
  #  assign("Knn", original_knn, envir = .GlobalEnv)
  #})

  # Test with invalid x structure
  invalid_x <- list(InvalidDF = processed$IntensityDF, CorrespondenceMatrix = processed$CorrespondenceMatrix, SpatialCoords = processed$SpatialCoords)
  expect_error(computeGearysC(invalid_x, verbose = FALSE, update_correspondence = FALSE))

  # Test with edge cases for coords
  processed_same_coords <- processed
  processed_same_coords$SpatialCoords <- data.frame(x = rep(1, nrow(processed_same_coords$SpatialCoords)), y = rep(1, nrow(processed_same_coords$SpatialCoords)))
  expect_silent(cur_test_same_coords <- computeGearysC(processed_same_coords, verbose = FALSE))

  processed_clustered_coords <- processed
  processed_clustered_coords$SpatialCoords <- data.frame(x = c(rep(1, nrow(processed_clustered_coords$SpatialCoords) / 2), rep(2, nrow(processed_clustered_coords$SpatialCoords) / 2)), y = rep(1, nrow(processed_clustered_coords$SpatialCoords)))
  expect_silent(cur_test_clustered_coords <- computeGearysC(processed_clustered_coords, verbose = FALSE))
})

