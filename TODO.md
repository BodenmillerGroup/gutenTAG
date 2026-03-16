- **make sure that untargeted metapeaks that are weirdly actually assigned to markers are removed**

- **Replace all mentions of "density" with "sparsity"**



---

## Code Review Findings

### Priority 1 — Bugs that will break on clean install or produce silent wrong results

- [ ] **Investigate `'x' should be a list object` error in `assignMetapeaks`**
  Users hitting `Error in .valid.assignMetapeaks(x, pre, refList, mz_threshold): 'x' should be a list object. It must explicitly be the output of the 'metapeakGeneration' function.` — determine whether this is a validity check that is too strict, a documentation issue, or a genuine misuse of the function.

- [ ] **Determine minimum compatible Cardinal version**
  Currently pinned to `>= 3.6.2` in DESCRIPTION to avoid a version mismatch (Cardinal 3.x changed `summarizeFeatures` API: `stat=` → `FUN=`, dropping `"nnzero"`). Investigate the oldest Cardinal version that works correctly with gutenTAG and update the version pin accordingly.
  **Note:** Cardinal 3.6.2 can be installed on R 4.3.1 (without upgrading to R 4.4) by installing `matter >= 2.6.2` and `Cardinal` directly via `remotes::install_github()`. Standard `BiocManager::install()` on Bioc 3.17/3.18 is insufficient as it ships older versions of these packages.

- [ ] **Move discarded duplicate peaks into `Untargeted` in `.removeDuplicates`**
  `utils.R`: peaks dropped during deduplication are currently lost. They should be appended to `x$Untargeted$UntargetedIntensity` and `x$Untargeted$UntargetedCorrespondence` so the output is consistent — targeted peaks that lost the duplicate contest become untargeted rather than disappearing entirely. Identify discarded peaks as those in `dup_correspondence$mz_location` not in `correct_peaks`, then cbind/rbind into the Untargeted fields.



- [ ] **Fix NaN propagation in `computeGearysC` for zero-variance channels**
  `computeGearysC.R` line 41: `var(X) * N` for a constant (all-zero) column produces 0, making `0/0 = NaN` in the correspondence matrix. Detect zero-variance channels before the formula and assign `NA` explicitly.


### Priority 2 — Convention violations and R CMD CHECK failures

- [ ] **Fix `MulticoreParam` Windows incompatibility in `preProcess.R`**
  `preProcess.R` line 30 unconditionally uses `MulticoreParam` (fork-based, fails on Windows). Apply the same conditional fallback pattern already used in `peakDetection.R`.


### Priority 3 — Performance and style issues

- [ ] **Replace vector-growing loops with pre-allocated vectors**
  Replace `c()` growing inside loops with `numeric(n)` + index assignment in:
  - `utils.R` lines 100–108 (`intravariance_vector`) and lines 184/193 (`bin_vector`)
  - `computeGearysC.R` lines 26/44 (`geary_vector`)
  Also replace `1:length(x)` / `1:ncol(df)` with `seq_along()` / `seq_len()` in `utils.R` line 185, `computeGearysC.R` line 27.

- [ ] **Replace `== TRUE` / `== T` comparisons throughout codebase**
  Use logical values directly. `T`/`F` abbreviations are especially dangerous. Affected: `utils.R` line 9, `readPanel.R` line 19, `computeGearysC.R` line 49, `asMSImagingExperiment.R` line 47, `assignMetapeaks.R`.

- [ ] **Add input validation to `generateSeedMz` for `density` parameter**
  Zero or negative `density` causes a cryptic error from `splus2R::peaks`. Add explicit validation and a `@param density` Roxygen entry.

### Priority 4 — Low-priority cleanup

- [ ] **Add QC plot for mean-variance residuals**
  Create a ggplot of `log(1 + mean_intensity)` vs `corrected_sd` residuals. Colour untargeted metapeaks in grey and annotated (targeted) ones in red. Useful for validating the `corrected_sd` metric in the correspondence matrix.

- [ ] **Remove dead and commented-out code blocks**
  - `smoothPeakCounts.R` lines 21–44: commented-out TODO implementations
  - `validityChecks.R` lines 188–382: commented-out validity functions (`.valid.pca`, `.valid.computeNMF`, etc.)
  - `utils.R` line 34: fix misleading comment and remove dead `endsWith(new_string, "-")` check on line 36
  - Add `*.swp` to `.gitignore`


- [ ] **Fix `imageChannel.R` `@return` documentation**
  Function plots as a side effect and returns the `cimg` object invisibly. Document this correctly.

### Runtime Bottlenecks — `assignMetapeaks.R`

- [ ] **Replace `cbind` growing loop in `.generateFinalIntensityDF`**
  `assignMetapeaks.R` lines 117–125: `cbind(final_intensity, intensity_temp)` copies the entire matrix on every iteration. Pre-allocate `matrix(0, nrow = n_pixels, ncol = nrow(correspondence))` and fill by column index. Likely the largest runtime bottleneck on big datasets.

- [ ] **Replace column loop in `.filterTIC` with index-based assignment**
  `assignMetapeaks.R` lines 239–248: `for` loop over columns using `!channel %in% filtered_vals` is O(n²) per column. Replace with direct index assignment using the already-computed `indices` vector. Also fixes the existing value-membership correctness bug (P1).

- [ ] **Vectorise `cor` computation in `.finaliseCorrespondence`**
  `assignMetapeaks.R` line 183: `apply(..., FUN = function(x) cor(...))` computes correlations column-by-column. Replace with a single `cor()` call on the full matrix.

