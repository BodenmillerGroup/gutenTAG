# TODO revisit
# TODO this tests the current behavior, but it would make sense to use non-circular convolution
test_that("smoothPeakCounts works when odd", {
  count_df <- data.frame(
    mz = c(1, 2, 3, 4, 5, 6, 7, 8, 9),
    count = c(10, 10, 10, 0, 0, 0, 10, 10, 10)
  )
  smooth_count_df <- smoothPeakCounts(count_df)
  expect_equal(
    smooth_count_df$count,
    c(17.5, 17.7, 17.7, 17.5, 13.9, 3.9, 0.4, 3.9, 13.9),
    tolerance = 0.1
  )
})


test_that("smoothPeakCounts works when even", {
  count_df <- data.frame(
    mz = c(1, 2, 3, 4, 5, 6, 7, 8),
    count = c(10, 10, 10, 0, 0, 0, 10, 10)
  )
  smooth_count_df <- smoothPeakCounts(count_df)
  expect_equal(smooth_count_df$count, c(16.6, 17.7, 17.7, 16.6, 8.9, 1.1, 1.1, 8.9), tolerance=0.1)
})
