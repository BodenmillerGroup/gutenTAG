test_that("generateSeedMz works", {
  count_df <- data.frame(
    mz = c(1, 2, 3, 4, 5, 6, 7, 8, 9),
    count = c(0, 200, 0, 2, 50, 80, 3, 0, 0)
  )
  seeds <- generateSeedMz(count_df)
  expect_equal(seeds, c(2, 6))
})
