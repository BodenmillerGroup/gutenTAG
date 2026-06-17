test_that("estimateMetapeaks: returns internally-consistent output on example data", {

  # use the same fixture the function's own @examples uses
  rdata_path <- system.file("extdata/Example_data/Example_processed.Rdata", package = "gutenTAG")
  load(rdata_path)

  count_df        <- results$metapeaks$count_df
  smooth_count_df <- results$metapeaks$count_smooth_df
  seed_mz         <- results$metapeaks$seed_mz

  res <- estimateMetapeaks(count_df, smooth_count_df, seed_mz = seed_mz)

  n <- max(res$propagation_selection)
  # every per-metapeak accumulator has one entry per metapeak, in order
  expect_length(res$metapeaks$center, n)
  expect_length(res$metapeaks$max, n)
  expect_length(res$metapeaks$width, n)
  expect_equal(dim(res$metapeaks$limits), c(n, 2L))
  # limits bracket the metapeak max
  expect_true(all(res$metapeaks$limits[, 1] <= res$metapeaks$max))
  expect_true(all(res$metapeaks$max <= res$metapeaks$limits[, 2]))

})

test_that("estimateMetapeaks: grouped reduction matches a hand reference", {

  # A single Gaussian bump with one seed segments into exactly one metapeak,
  # so the grouped-reduction output must equal the plain weighted mean / argmax
  # over the whole histogram.
  mz    <- seq(100, 110, by = 0.1)
  count <- exp(-((mz - 105)^2) / (2 * 0.5^2)) * 100
  count_df <- data.frame(mz = mz, count = count)

  res <- estimateMetapeaks(count_df = count_df, smooth_count_df = count_df,
                           seed_mz = 105, detection_threshold = 0)

  expect_equal(max(res$propagation_selection), 1)
  expect_equal(res$metapeaks$center,
               stats::weighted.mean(mz, count), tolerance = 1e-10)
  expect_equal(res$metapeaks$max, mz[which.max(count)], tolerance = 0)

})

test_that("estimateMetapeaks: single-metapeak case yields correctly-shaped output", {

  # Guards the cbind/tapply reassembly: with one metapeak, limits must still be
  # a 1x2 matrix (not a dropped numeric vector) and the scalar accumulators
  # must be length 1.
  mz    <- seq(100, 110, by = 0.1)
  count <- exp(-((mz - 105)^2) / (2 * 0.5^2)) * 100
  count_df <- data.frame(mz = mz, count = count)

  res <- estimateMetapeaks(count_df = count_df, smooth_count_df = count_df,
                           seed_mz = 105, detection_threshold = 0)

  # exactly one metapeak
  expect_equal(max(res$propagation_selection), 1)

  # limits reassembled as a 1x2 matrix
  expect_true(is.matrix(res$metapeaks$limits))
  expect_equal(dim(res$metapeaks$limits), c(1L, 2L))

  # scalar accumulators are length 1
  expect_length(res$metapeaks$center, 1)
  expect_length(res$metapeaks$max, 1)
  expect_length(res$metapeaks$width, 1)

  # the single metapeak is located on the bump
  expect_equal(res$metapeaks$max, 105, tolerance = 0.05)

})

test_that("estimateMetapeaks: fixed.limits sets a symmetric window around the max", {

  mz    <- seq(100, 110, by = 0.1)
  count <- exp(-((mz - 105)^2) / (2 * 0.5^2)) * 100
  count_df <- data.frame(mz = mz, count = count)

  res <- estimateMetapeaks(count_df = count_df, smooth_count_df = count_df,
                           seed_mz = 105, detection_threshold = 0,
                           fixed.limits = 0.5)

  expect_equal(res$metapeaks$limits[1, 1], res$metapeaks$max - 0.5, tolerance = 0)
  expect_equal(res$metapeaks$limits[1, 2], res$metapeaks$max + 0.5, tolerance = 0)

})
