- **make sure that untargeted metapeaks that are weirdly actually assigned to markers are removed**

- **Replace all mentions of "density" with "sparsity"**



---

## Code Review Findings

### Priority 1 — Bugs that will break on clean install or produce silent wrong results


- [ ] **Determine minimum compatible Cardinal version**
  Currently pinned to `>= 3.6.2` in DESCRIPTION to avoid a version mismatch (Cardinal 3.x changed `summarizeFeatures` API: `stat=` → `FUN=`, dropping `"nnzero"`). Investigate the oldest Cardinal version that works correctly with gutenTAG and update the version pin accordingly.
  **Note:** Cardinal 3.6.2 can be installed on R 4.3.1 (without upgrading to R 4.4) by installing `matter >= 2.6.2` and `Cardinal` directly via `remotes::install_github()`. Standard `BiocManager::install()` on Bioc 3.17/3.18 is insufficient as it ships older versions of these packages.


- [ ] **Move discarded duplicate peaks into `Untargeted` in `.removeDuplicates`**
  `utils.R`: peaks dropped during deduplication are currently lost. They should be appended to `x$Untargeted$UntargetedIntensity` and `x$Untargeted$UntargetedCorrespondence` so the output is consistent — targeted peaks that lost the duplicate contest become untargeted rather than disappearing entirely. Identify discarded peaks as those in `dup_correspondence$mz_location` not in `correct_peaks`, then cbind/rbind into the Untargeted fields.





### Priority 2 — Convention violations and R CMD CHECK failures

- [ ] **Fix `MulticoreParam` Windows incompatibility in `preProcess.R`**
  `preProcess.R` line 30 unconditionally uses `MulticoreParam` (fork-based, fails on Windows). Apply the same conditional fallback pattern already used in `peakDetection.R`.


### Priority 3 — Performance and style issues


### Priority 4 — Low-priority cleanup


- [ ] **Add QC plot for mean-variance residuals**
  Create a ggplot of `log(1 + mean_intensity)` vs `corrected_sd` residuals. Colour untargeted metapeaks in grey and annotated (targeted) ones in red. Useful for validating the `corrected_sd` metric in the correspondence matrix.

- [ ] **Remove dead and commented-out code blocks**
  - `smoothPeakCounts.R` lines 21–44: commented-out TODO implementations
  - `validityChecks.R` lines 188–382: commented-out validity functions (`.valid.pca`, `.valid.computeNMF`, etc.)
  - `utils.R` line 34: fix misleading comment and remove dead `endsWith(new_string, "-")` check on line 36
  - Add `*.swp` to `.gitignore`




