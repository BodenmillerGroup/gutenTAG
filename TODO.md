- **make sure that untargeted metapeaks that are weirdly actually assigned to markers are removed**

---

## Code Review Findings

### Priority 1 — Bugs that will break on clean install or produce silent wrong results


- [ ] **Determine minimum compatible Cardinal version**
  Currently pinned to `>= 3.6.2` in DESCRIPTION to avoid a version mismatch (Cardinal 3.x changed `summarizeFeatures` API: `stat=` → `FUN=`, dropping `"nnzero"`). Investigate the oldest Cardinal version that works correctly with gutenTAG and update the version pin accordingly.
  **Note:** Cardinal 3.6.2 can be installed on R 4.3.1 (without upgrading to R 4.4) by installing `matter >= 2.6.2` and `Cardinal` directly via `remotes::install_github()`. Standard `BiocManager::install()` on Bioc 3.17/3.18 is insufficient as it ships older versions of these packages.


- [ ] **Move discarded duplicate peaks into `Untargeted` in `.removeDuplicates`**
  `utils.R`: peaks dropped during deduplication are currently lost. They should be appended to `x$Untargeted$UntargetedIntensity` and `x$Untargeted$UntargetedCorrespondence` so the output is consistent — targeted peaks that lost the duplicate contest become untargeted rather than disappearing entirely. Identify discarded peaks as those in `dup_correspondence$mz_location` not in `correct_peaks`, then cbind/rbind into the Untargeted fields.
  **Note:** `plotMeanVarianceResiduals` (`R/plotQC.R`) reads directly from `x$Untargeted$UntargetedIntensity` and matches column names against `CorrespondenceMatrix$mz_location`. Any changes to how untargeted metapeaks are stored (structure, naming, or which peaks are included) will require updating this function accordingly.





### Priority 2 — Convention violations and R CMD CHECK failures


### Priority 3 — Performance and style issues

- [ ] **Add BiocParallel parallelisation to `generateMetapeaks` and `assignMetapeaks`**
  Plan saved at `~/.claude/plans/generic-bouncing-scone.md`. Parallelise the `for` loop over metapeaks in `estimateMetapeaks` (lines 37–76) and the `for` loop in `.generateFinalIntensityDF` (lines 136–141) using `BiocParallel::bplapply`. Add `BPPARAM = BiocParallel::bpparam()` to both exported functions and their internal helpers. Add validity checks and `@param BPPARAM` roxygen entries. Add serial-vs-parallel equality tests using `SnowParam(2)`.


### Priority 4 — Low-priority cleanup


- [ ] **Remove dead and commented-out code blocks**
  - `smoothPeakCounts.R` lines 21–44: commented-out TODO implementations
  - `validityChecks.R` lines 188–382: commented-out validity functions (`.valid.pca`, `.valid.computeNMF`, etc.)
  - `utils.R` line 34: fix misleading comment and remove dead `endsWith(new_string, "-")` check on line 36
  - Add `*.swp` to `.gitignore`




