- **make sure that untargeted metapeaks that are weirdly actually assigned to markers are removed**



---

## Code Review Findings

### Priority 1 — Bugs that will break on clean install or produce silent wrong results

- [ ] **Fix `.removeDuplicates` Inf/-Inf crash on empty shifted mz vector**
  `utils.R` lines 152–162: when `shift == "right"`, `correct_shifted_mzs` can be `numeric(0)`, making `min(numeric(0))` return `Inf` — silently dropping that marker. Add a guard before calling `min()`/`max()` on a possibly empty vector.

- [ ] **Fix `.filterTIC` using value-membership instead of index-membership**
  `assignMetapeaks.R` lines 243–246: `!channel %in% filtered_vals` incorrectly retains pixels that share intensity values with passing pixels. Fix: use index-based filtering (`channel[-passing_indices] <- NA`).

- [ ] **Fix `corrected_sd` computation in `.finaliseCorrespondence`**
  `assignMetapeaks.R` line 173: `residuals + sd_intensity / sd_intensity` simplifies to `residuals + 1` (sd/sd == 1, NaN when sd == 0). Likely unintentional — investigate the intended normalization and fix.

- [ ] **Fix NaN propagation in `computeGearysC` for zero-variance channels**
  `computeGearysC.R` line 41: `var(X) * N` for a constant (all-zero) column produces 0, making `0/0 = NaN` in the correspondence matrix. Detect zero-variance channels before the formula and assign `NA` explicitly.

- [ ] **Fix `.otsu_thresholding` NA propagation from `var()` on short vectors**
  `utils.R` lines 105–106: `var()` returns `NA` for length-0 or length-1 vectors, causing `NA` in `intravariance_vector` and `which.min()` returning the wrong bin. Fix with `if (length(x) < 2) NA else var(x)`.

### Priority 2 — Convention violations and R CMD CHECK failures

- [ ] **Fix discarded `relocate()` call in `.cleanPanel`**
  `utils.R` line 14: `dplyr::relocate(panel, ...)` return value is not assigned back to `panel`, making it a dead call that can produce wrong column ordering. Remove or fix the assignment.

- [ ] **Fix `MulticoreParam` Windows incompatibility in `preProcess.R`**
  `preProcess.R` line 30 unconditionally uses `MulticoreParam` (fork-based, fails on Windows). Apply the same conditional fallback pattern already used in `peakDetection.R`.

- [ ] **Fix dead code and typo in `.valid.peakDetection`**
  `validityChecks.R` lines 58–62: unreachable `isS4(x)` branch and typo `"MSImagingExperiement"` in the error message.

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

- [ ] **Remove dead and commented-out code blocks**
  - `smoothPeakCounts.R` lines 21–44: commented-out TODO implementations
  - `validityChecks.R` lines 188–382: commented-out validity functions (`.valid.pca`, `.valid.computeNMF`, etc.)
  - `utils.R` line 34: fix misleading comment and remove dead `endsWith(new_string, "-")` check on line 36
  - Add `*.swp` to `.gitignore`

- [ ] **Replace placeholder test in `test-plotPixels.R`**
  File contains only `expect_equal(2 * 2, 4)`. Write real tests for `plotPixels` or remove the file.

- [ ] **Fix `imageChannel.R` `@return` documentation**
  Function plots as a side effect and returns the `cimg` object invisibly. Document this correctly.

