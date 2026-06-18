# Expensive fixture for the large-n Geary's C test.
#
# Bioconductor rule: tests must not hit the network and expensive setup belongs in
# a helper-*.R file rather than inline in a test. This builds a ~72,000-pixel
# regular grid that CROSSES the 2^16 (65,536) pixel boundary at which the old
# dense as.matrix(dist(coords)) builder overflowed R's 32-bit vector index.
#
# The returned object is a minimal "processed" list — the only structure
# computeGearysC() and its validity check require (IntensityDF, SpatialCoords,
# CorrespondenceMatrix) — not a full assignMetapeaks() output.

make_largen_processed <- function(nx = 360, ny = 200) {

  # column-major lattice on a regular integer grid starting at 1
  coords <- expand.grid(x = seq_len(nx), y = seq_len(ny))

  set.seed(1)

  # channel 1: spatially structured (smooth gradient in x + small noise) ->
  #   strong positive spatial autocorrelation -> LOW Geary's C
  # channel 2: pure spatial noise -> no structure -> Geary's C near 1
  intensity <- data.frame(
    structured = coords$x + rnorm(nrow(coords), sd = 0.01),
    noise      = rnorm(nrow(coords))
  )

  correspondence <- data.frame(
    marker      = colnames(intensity),
    mz_location = c(900.1, 901.2)
  )

  list(
    IntensityDF          = intensity,
    SpatialCoords        = coords,
    CorrespondenceMatrix = correspondence
  )
}
