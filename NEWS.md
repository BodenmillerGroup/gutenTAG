# gutenTAG 0.99.17

* Added `plotMassShiftVsIntensity()`, a QC scatter of each targeted peak's
  signed mass shift (`mz_location - expected_mz_location`) against its mean
  intensity on a log10 axis, as a screen for falsely-assigned markers. Points
  whose absolute shift exceeds `shift_threshold` (defaulting to `2 * sd` of the
  shifts) are flagged in a contrasting colour and labelled with their marker
  name via `ggrepel`. Customisable via `points`, `labels`, `shift_threshold`,
  `show_labels` and `palette`.
* `plotIntensityDistribution()` gained a `transformation` argument
  (`"none"`, `"log"`, `"z-scaled"`) for viewing per-marker intensity
  distributions on the raw scale, the log10 scale (with a `+1` pseudocount so
  zero-intensity pixels are handled), or as the z-score of the log10 intensity
  pooled across all markers. Marker order on the x-axis is fixed by raw mean
  intensity so the three views are directly comparable.
* `computeGearysC()` now scales to large images. Its internal spatial-weight
  matrix builder was rewritten to construct the rook (`k = 4`) adjacency
  directly from the integer pixel lattice in O(n) time and memory, instead of
  materialising a dense n x n pairwise-distance matrix that overflowed R's
  vector index above ~65,000 pixels. Note a deliberate, minor behaviour change:
  edge and hole pixels now connect only to their existing 2-3 rook neighbours
  rather than 4 force-picked nearest points (which previously pulled in
  diagonals at boundaries). This is the geometrically correct neighbourhood and
  slightly shifts edge-pixel scores.

# gutenTAG 0.99.15

* Initial Bioconductor submission.
* Added `plotMassShift()`, a QC plot of the deviation between expected and
  picked m/z, with `"distribution"` (beeswarm) and `"spectrum"` (along the m/z
  axis) views.
* `plotMetapeaks()` now draws a dashed line at the detection threshold; the
  absolute threshold is stored in `generateMetapeaks()` output as
  `params$detection_threshold`.
* `plotMeanVariance()` now draws a linear regression line (customisable via the
  `regression` argument, e.g. `regression = list(color = "blue", se = TRUE)`)
  and annotates its R-squared and slope in the top-left corner (toggle with
  `show_stats`, style with `stats`). The line and stats are fitted over the
  unannotated points, or over the targeted points when `annotated_only = TRUE`.
  The smoothing method is configurable (`regression = list(method = "loess")`,
  etc.; `"gam"` and `"rlm"` require the suggested `mgcv`/`MASS` packages); the
  stats box refits and labels that same method, reporting slope only when the
  fit is linear in the predictor.
* `plotMeanVsSNR()` gained the same configurable regression line and
  method-coupled stats box (`regression`, `show_stats`, `stats` arguments),
  fitted on the linear mean/SNR axes.
* All QC plot functions except `plotMetapeaks()` now accept per-layer aesthetic
  lists for flexible styling of individual plot layers (e.g. `points`,
  `unannotated_points`, `bar`, `violin`, `boxplot`, `histogram`, `line`,
  `vline`, `rect`, `labels`), e.g. `points = list(size = 5, fill = "red")`,
  `line = list(color = "blue")`. Plots whose colour/fill encodes a data grouping
  (`plotSNR()`, `plotSNRHistogram()`, `plotSpatialSNR()`, `plotQCOverview()`,
  `plotTICHistogram()`, `plotTICSpatial()`) additionally accept a `palette`
  argument to recolour the groups, e.g. `palette = c(Signal = "purple")`.
* Exported `removeDuplicates()` (previously the internal `.removeDuplicates`),
  with documentation.
* Excluded rendered QC report artifacts (`inst/qc_report_files/`, `*.html`) from
  the build via `.Rbuildignore`, reducing installed package size by ~17 MB.
