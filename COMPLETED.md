# Completed Tasks

## Version 0.12.1 — Branch: bioc-fixes

- make sure default parameter for watershed threshold is 0.01
- make sure processed$IntensityDF is ordered by m/z value
- make sure that CorrespondenceMatrix$expected_mz_location is numeric and not character
- untargeted metapeak columns should be assigned the name of that m/z value (as.character)
- Fix EBImage missing from DESCRIPTION and NAMESPACE — already present: EBImage in DESCRIPTION Imports, importFrom(EBImage,Image) and importFrom(EBImage,propagate) in NAMESPACE, and @importFrom tags in asCytoImageList.R and segmentPeakCounts.R
- Fix abind missing from DESCRIPTION and NAMESPACE — already present: abind in DESCRIPTION Imports, importFrom(abind,abind) in NAMESPACE, and @importFrom abind abind in asCytoImageList.R
- Fix DESCRIPTION: Author format, License, LazyData, BiocStyle placement — already correct; removed erroneous @import BiocStyle from preProcess.R and import(BiocStyle) from NAMESPACE
- Fix duplicate Cardinal import in NAMESPACE — no duplicate present; @import Cardinal was already absent from generateMetapeaks.R
- README.md overhaul: rewrote opening paragraph; fixed typos; fixed heading hierarchy (H4→H3); standardised code fences; fixed preProcess description (smoothSpectra→smooth); fixed assignMetapeaks argument order; fixed assignMetapeaks output list (added SummarySpectra, named Untargeted sub-elements); updated .ibd file note; resolved metapeakGeneration→generateMetapeaks; added QC and conversion sections; added Citation, License, Contributing sections; added Bioconductor install instructions with BiocManager note; added status badge placeholders; removed PCA/NMF/UMAP direct-function claim; expanded cytoviewer aside
- Fix .valid.generateMetapeaks to validate x is MSImagingExperiment
- Fix dead code and typo in .valid.peakDetection — removed unreachable isS4(x) branch and misspelled error message from validityChecks.R
- Remove placeholder test-plotPixels.R — file contained only usethis boilerplate; plotPixels is a cytomapper function, not a gutenTAG export
- Fix .otsu_thresholding NA propagation — replaced var() with safe_var() helper that returns 0 for partitions with fewer than 2 elements, preventing NA in intravariance_vector and wrong which.min() result
- Fix .removeDuplicates: refactored to accept x directly instead of sample$processed; pre-allocate correct_peaks; add shift input validation; guard against empty correct_shifted_mzs with fallback and warning; collapse symmetric left/right branches
- Fix .filterTIC value-membership bug — replaced !channel %in% filtered_vals with channel[-indices] <- NA for correct index-based pixel filtering; note this is a behaviour change for datasets where failing pixels share intensity values with passing pixels
- Fix corrected_sd computation in .finaliseCorrespondence — changed from residuals + sd_intensity / sd_intensity (which simplified to residuals + 1) to residuals / sd_intensity (dimensionless, comparable across markers)
- Update tests for imageChannel — replaced imager-based tests (cimg class, old parameter names, hard-coded pixel values) with ggplot2-based tests covering return type, list and df+coords input, channel by index/name, quantile thresholding, palette, na_colour, validity checks, and theme passthrough
- Replace `== TRUE` / `== T` comparisons — removed redundant `== TRUE` in `computeGearysC.R:50` and `asMSImagingExperiment.R:46`; no `== T` instances found
- Fix NaN propagation in `computeGearysC` for zero-variance channels — added `is.na(Var_X) || Var_X == 0` guard returning `NA_real_`; switched to `var(X, na.rm = TRUE)` to handle pixel-level NAs silently; replaced `1:ncol(df)` with `seq_len(ncol(df))`; updated roxygen docs with full `@param`, corrected `@return` for both branches, and added `@examples`; fixed misleading validity error message
- Replace vector-growing loops with pre-allocated vectors — `geary_vector <- c()` replaced with `numeric(ncol(df))` and index assignment in `computeGearysC.R`; `intravariance_vector` and `bin_vector` in `utils.R` were already pre-allocated with `numeric()`/`logical()` and `seq_len`; strengthened `.valid.computeGearysC` to check for required list fields
- Add `shift = "closest"` option to `.removeDuplicates` — added `.choose_closest()` helper; updated input guard to accept `"closest"`; inserted `"closest"` branch calling `.choose_closest()` directly (no directional filter, no fallback warning); added `test_that` blocks for `.choose_closest` (5 cases) and `.removeDuplicates` (6 cases) in `test-utils.R`
