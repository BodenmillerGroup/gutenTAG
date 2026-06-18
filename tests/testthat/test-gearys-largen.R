# Large-n regression test for the grid-native rook adjacency rewrite.
#
# The fixture (see helper-gearys-largen.R) builds ~72,000 pixels, crossing the
# 65,536-pixel boundary at which the OLD dense as.matrix(dist(coords)) builder
# overflowed R's 32-bit vector index and errored / OOMed. These tests confirm the
# new O(n) rook builder runs to completion and produces sensible scores.
#
# skip_on_cran(): the fixture is deliberately large and slow to build.

test_that("computeGearysC scales past the 2^16 pixel overflow boundary", {
  skip_on_cran()

  big <- make_largen_processed(nx = 360, ny = 200)

  # sanity: we really did cross the boundary that broke the old code
  expect_gt(nrow(big$SpatialCoords), 70000)

  # the OLD as.matrix(dist(coords)) code errored / OOMed here; the new builder
  # must run without error.
  expect_no_error(res <- computeGearysC(big, verbose = FALSE))

  # one score per marker channel, all finite
  expect_equal(length(res), 2)
  expect_true(all(is.finite(res)))

  # the spatially structured channel must score BELOW the pure-noise channel
  # (lower Geary's C == stronger positive spatial autocorrelation)
  expect_lt(res[1], res[2])
})

test_that(".knn_weight_matrix stays sparse at large n (no re-densification)", {
  skip_on_cran()

  big <- make_largen_processed(nx = 360, ny = 200)
  coords <- as.matrix(big$SpatialCoords)
  n <- nrow(coords)

  W <- gutenTAG:::.knn_weight_matrix(coords, k = 4)

  expect_s4_class(W, "sparseMatrix")

  # rook adjacency stores at most 4 edges per pixel (< 20n is a loose guard).
  # This assertion fails LOUDLY if anyone reintroduces a dense n x n build.
  expect_lt(Matrix::nnzero(W), 20 * n)
})
