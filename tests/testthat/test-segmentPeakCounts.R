test_that("segmentPeakCounts works", {
  count_df <- data.frame(
    mz = c(10, 20, 30, 40, 50, 60, 70, 80, 90),
    count = c(20, 19, 20, 0, 0, 19, 20, 19, 0)
  )
  seed_mz <- c(20, 60)
  segmentation <- segmentPeakCounts(count_df, seed_mz, detection_threshold = 0)
  expect_equal(segmentation, c(1, 1, 1, 0, 0, 2, 2, 2, 0))
})
