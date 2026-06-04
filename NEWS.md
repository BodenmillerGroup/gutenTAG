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
