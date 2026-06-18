# Tests for the internal rook-adjacency spatial weight builder .knn_weight_matrix.
#
# This helper was rewritten from a dense force-k-nearest implementation to a
# grid-native rook adjacency. The structural invariants below (symmetry, zero
# diagonal, unit weights, exact neighbour counts) are the contract of the new
# implementation. We do NOT pin the old byte-for-byte squared-distance values —
# those intentionally changed (see the corner/edge count assertions, which
# document the sanctioned edge-behaviour change vs. the old force-4-nearest code).

test_that(".knn_weight_matrix builds a symmetric unit-weight rook graph (3x3)", {

  # expand.grid is column-major: x cycles fastest, so pixel indices are
  #   1:(1,1) 2:(2,1) 3:(3,1) 4:(1,2) 5:(2,2) 6:(3,2) 7:(1,3) 8:(2,3) 9:(3,3)
  coords <- as.matrix(expand.grid(x = 1:3, y = 1:3))

  W <- gutenTAG:::.knn_weight_matrix(coords, k = 4)

  # structural invariants
  expect_s4_class(W, "sparseMatrix")
  expect_true(Matrix::isSymmetric(W))
  expect_equal(Matrix::diag(W), rep(0, nrow(coords)))

  # unit weights: every stored edge carries value 1
  nz <- W@x
  expect_true(all(nz == 1))

  # --- EXACT rook adjacency (the intended change vs. old force-4-nearest) ---
  # neighbour count per pixel == number of stored entries in its row
  deg <- Matrix::rowSums(W != 0)

  # centre pixel (2,2) is index 5 — exactly 4 rook neighbours
  expect_equal(deg[5], 4)

  # corner pixels (1,1),(3,1),(1,3),(3,3) are indices 1,3,7,9 — exactly 2 each.
  # The OLD code would have forced 4 neighbours here, pulling in diagonals; the
  # rook graph correctly connects a corner to only its 2 existing neighbours.
  expect_equal(unname(deg[c(1, 3, 7, 9)]), c(2, 2, 2, 2))

  # edge (non-corner) pixels (2,1),(1,2),(3,2),(2,3) are indices 2,4,6,8 — 3 each
  expect_equal(unname(deg[c(2, 4, 6, 8)]), c(3, 3, 3, 3))

  # spot-check the centre's actual neighbours: (2,2)=5 links to
  # (1,2)=4, (3,2)=6, (2,1)=2, (2,3)=8
  expect_equal(sort(which(W[5, ] != 0)), c(2, 4, 6, 8))
})

test_that(".knn_weight_matrix invariants hold on a 5x5 grid", {

  coords <- as.matrix(expand.grid(x = 1:5, y = 1:5))
  n <- nrow(coords)

  W <- gutenTAG:::.knn_weight_matrix(coords, k = 4)

  expect_s4_class(W, "sparseMatrix")
  expect_true(Matrix::isSymmetric(W))
  expect_equal(Matrix::diag(W), rep(0, n))
  expect_true(all(W@x == 1))

  deg <- Matrix::rowSums(W != 0)

  # interior pixel (3,3): column-major index = 3 + (3-1)*5 = 13 -> 4 neighbours
  expect_equal(deg[13], 4)

  # four corners: indices 1,5,21,25 -> 2 neighbours each
  expect_equal(unname(deg[c(1, 5, 21, 25)]), c(2, 2, 2, 2))

  # a 5x5 grid has 4 corners (deg 2), 12 border non-corners (deg 3), 9 interior
  # (deg 4). Check the degree distribution matches.
  expect_equal(sum(deg == 2), 4)
  expect_equal(sum(deg == 3), 12)
  expect_equal(sum(deg == 4), 9)

  # total stored entries == sum of degrees == 2 * (number of undirected edges).
  # A 5x5 grid has 5*4 horizontal + 4*5 vertical = 40 undirected edges -> 80 nnz.
  expect_equal(Matrix::nnzero(W), 80)
})

test_that(".knn_weight_matrix connects only existing neighbours across a hole", {

  # a 3x3 grid with the centre pixel (2,2) removed — simulates a tissue hole.
  full   <- expand.grid(x = 1:3, y = 1:3)
  coords <- as.matrix(full[!(full$x == 2 & full$y == 2), ])
  n <- nrow(coords)               # 8 pixels

  W <- gutenTAG:::.knn_weight_matrix(coords, k = 4)

  expect_true(Matrix::isSymmetric(W))
  expect_equal(Matrix::diag(W), rep(0, n))
  expect_true(all(W@x == 1))

  # the pixel (2,1) used to have (2,2) as a neighbour; with the hole it keeps only
  # (1,1) and (3,1) — i.e. 2 neighbours, NOT a forced 4. This is the sanctioned
  # change: holes/edges no longer pull in non-adjacent pixels.
  row_21 <- which(coords[, "x"] == 2 & coords[, "y"] == 1)
  expect_equal(sum(W[row_21, ] != 0), 2)
})
