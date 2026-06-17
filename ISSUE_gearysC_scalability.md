# Issue: `computeGearysC` does not scale — OOM and integer overflow on large samples

**Status:** open / unaddressed (documented 2026-06-17)
**Severity:** blocking for QC on real MALDI samples (≳65k pixels)
**Affected code:** `computeGearysC()` (`R/computeGearysC.R`) via the internal helper `.knn_weight_matrix()` (`R/utils.R:77`)

## Summary

`computeGearysC()` fails on real-sized MALDI samples. The spatial weight matrix it
depends on is built by materialising the **full dense `n × n` pairwise-distance
matrix** (`n` = number of pixels). This is O(n²) in both memory and time, and on
large images it fails in two different ways depending on the memory limit:

- **Out of memory** (smaller `--mem`): the dense matrix can't be allocated.
- **Integer overflow** (larger `--mem`): `dist()` produces `n(n−1)/2` values; once
  that exceeds `.Machine$integer.max` (~2.1e9, i.e. around **n ≈ 65,000 pixels**),
  the indexing in `as.matrix.dist()` overflows and the call fails.

Because the second failure is an overflow rather than a memory shortage, **adding
RAM cannot fix it** — the construction is fundamentally wrong for large `n`.

## Where it happens

`computeGearysC()` (`R/computeGearysC.R:54`) calls:

```r
spatial_weight_matrix <- .knn_weight_matrix(as.matrix(coords), k = 4)
```

`.knn_weight_matrix()` (`R/utils.R:77`) builds several dense `n × n` objects before
returning a sparse result:

```r
.knn_weight_matrix <- function(coords, k) {
  n  <- nrow(coords)
  d2 <- as.matrix(dist(coords))^2   # <-- dense n x n  (the problem)
  diag(d2) <- Inf
  knn_adj <- matrix(FALSE, n, n)    # <-- another dense n x n
  for (i in seq_len(n)) {
    knn_adj[i, order(d2[i, ])[seq_len(k)]] <- TRUE
  }
  knn_sym <- knn_adj | t(knn_adj)   # <-- and more n x n
  ...
  Matrix::sparseMatrix(...)         # only the OUTPUT is sparse
}
```

The output is a sparse k=4 ("rook") neighbour graph, but the intermediate
`as.matrix(dist(coords))` and `matrix(FALSE, n, n)` are dense.

## Memory cost (one dense `n × n` double matrix)

| pixels (n) | one n×n matrix | with the ~4 copies built | dist() length n(n−1)/2 |
|-----------:|---------------:|-------------------------:|-----------------------:|
| 30,000     | ~7 GB          | ~30–40 GB                | 4.5e8 (ok)             |
| 50,000     | ~20 GB         | ~80–100 GB               | 1.25e9 (ok)            |
| 65,000     | ~34 GB         | ~140 GB                  | ~2.1e9 (**overflow boundary**) |
| 100,000    | ~80 GB         | >250 GB                  | 5e9 (**overflows**)    |

## Observed failures (hyperplex_paper QC, sample B20_30928_A)

Rendering the QC report (which calls `computeGearysC(processed, update_correspondence = TRUE)`)
died at the `gearys_c_barchart` chunk:

- `--mem=128G`: SLURM `Detected 1 oom_kill event ... OOM Killed`.
- `--mem=256G`:
  ```
  Error in `sequence.default()`: ! 'from' contains NAs
   1. └─gutenTAG::computeGearysC(processed, update_correspondence = TRUE)
   2.   └─gutenTAG:::.knn_weight_matrix(as.matrix(coords), k = 4)
   3.     ├─base::as.matrix(dist(coords))
   ...
   6.       └─base::sequence.default(...)
  ```
  i.e. `as.matrix(dist(coords))` overflowed for this pixel count.

## Why it shouldn't need this

A k-nearest-neighbour graph needs only O(n·k) storage. For MALDI imaging the pixels
sit on a regular lattice and the `k = 4` "rook" neighbourhood is just the up/down/
left/right grid neighbours — buildable directly in O(n) without any pairwise-distance
matrix.

## Candidate fixes (for the plan; not yet implemented)

1. **Grid-native rook adjacency (no new dependency).** Build the sparse adjacency by
   coordinate lookup of `(x±1, y)` / `(x, y±1)` offsets. O(n) memory/time, exact for
   lattice data, base R + Matrix only. Caveat: at tissue edges/holes it connects a
   pixel to only its *existing* rook neighbours (2–3), whereas the current code
   force-picks 4 nearest including diagonals — a small, arguably-more-correct
   semantic change.

2. **Proper KNN search (new dependency).** Replace `as.matrix(dist())` with
   `FNN::get.knn(coords, k)` (kd-tree, O(n log n) time, O(n·k) memory). Preserves the
   exact "k nearest by Euclidean distance" semantics for arbitrary (non-grid) coords,
   but adds an `FNN` import that must be available in the install environment.

## Notes for whoever picks this up

- Geary's C is invariant to a global scaling of all weights, so using unit weights
  vs. squared-distance weights for a uniform-spacing grid gives the same score — the
  fix doesn't need to preserve the exact `d²` weight values, only the adjacency.
- Add a regression test with a known small grid (e.g. a 3×3 or 5×5 lattice) where the
  rook adjacency and resulting Geary's C can be computed by hand, plus a large-`n`
  smoke test (e.g. n > 70,000) that would have caught the overflow.
- Verify against `computeSNR` / the QC report path: only Geary's C builds spatial
  matrices, so this is the sole scalability blocker in the QC pipeline.
