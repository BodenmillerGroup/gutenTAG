- **make sure that untargeted metapeaks that are weirdly actually assigned to markers are removed**

- estimation criterion for what watershed threshold to use



---

## Code Review Findings

### Priority 1 — Bugs that will break on clean install or produce silent wrong results

- [ ] **Fix EBImage missing from DESCRIPTION and NAMESPACE**
  `EBImage::propagate` (segmentPeakCounts.R) and `EBImage::Image` (asCytoImageList.R) are called but EBImage is not in `DESCRIPTION` Imports or `NAMESPACE`. Add to `Imports` and add `@importFrom EBImage propagate Image` Roxygen tags.

- [ ] **Fix abind missing from DESCRIPTION and NAMESPACE**
  `abind()` called unqualified in asCytoImageList.R but not declared as a dependency. Add to `Imports` and add `@importFrom abind abind`.

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

- [ ] **Fix DESCRIPTION: Author format, License, LazyData, BiocStyle placement**
  - `Author` must use `Authors@R:` with `person()` objects
  - `License: MIT License` → `MIT + file LICENSE` (valid SPDX)
  - Remove `LazyData: true` (no `data/` directory exists)
  - Move `BiocStyle` from `Imports` to `Suggests`

- [ ] **Fix duplicate Cardinal import in NAMESPACE**
  Both `import(Cardinal, except = peaks)` and `import(Cardinal)` appear in NAMESPACE. Remove the unconditional `@import Cardinal` from `generateMetapeaks.R` to avoid R CMD CHECK warnings.

- [ ] **Fix dead code and typo in `.valid.peakDetection`**
  `validityChecks.R` lines 58–62: unreachable `isS4(x)` branch and typo `"MSImagingExperiement"` in the error message.

- [ ] **Fix `louvainCluster` silently overwriting `coords` parameter**
  `louvainCluster.R` line 51: validated `coords` parameter is unconditionally replaced by `x$SpatialCoords`. When `x` is a bare data frame, this yields `NULL`. Either use the validated parameter or remove it from the signature.

### Priority 3 — Performance and style issues

- [ ] **Replace vector-growing loops with pre-allocated vectors**
  Replace `c()` growing inside loops with `numeric(n)` + index assignment in:
  - `utils.R` lines 100–108 (`intravariance_vector`) and lines 184/193 (`bin_vector`)
  - `computeGearysC.R` lines 26/44 (`geary_vector`)
  Also replace `1:length(x)` / `1:ncol(df)` with `seq_along()` / `seq_len()` in `utils.R` line 185, `computeVariogram.R` line 23, `computeGearysC.R` line 27.

- [ ] **Replace `== TRUE` / `== T` comparisons throughout codebase**
  Use logical values directly. `T`/`F` abbreviations are especially dangerous. Affected: `utils.R` line 9, `readPanel.R` line 19, `computeGearysC.R` line 49, `asMSImagingExperiment.R` line 47, `assignMetapeaks.R`.

- [ ] **Add input validation to `generateSeedMz` for `density` parameter**
  Zero or negative `density` causes a cryptic error from `splus2R::peaks`. Add explicit validation and a `@param density` Roxygen entry.

- [ ] **Fix `.valid.generateMetapeaks` to validate `x` is `MSImagingExperiment`**
  `validityChecks.R` lines 67–101: `generateMetapeaks("test")` fails with a cryptic downstream error. Add an explicit `is(x, "MSImagingExperiment")` check.

- [ ] **Fix error messages in `.valid.computeVariogram` using wrong variable names**
  `validityChecks.R` lines 269 and 275: both messages say `"'x' should be a dataframe"` but the arguments are named `df` and `coords`.

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

---

## README.md Updates

> General note: when next editing README.md, strip out outdated or redundant content and tighten the prose throughout. The current README is roughly half of a complete workflow and contains several inaccuracies.

### Remove / correct outdated content

- [ ] **Remove or update the PCA/NMF/UMAP claim** in the introduction — no exported wrapper functions exist for these yet. Either add the functions or replace the claim with a note that conversion to `SpatialExperiment`/`AnnData` enables downstream dimensionality reduction.
- [ ] **Fix the `assignMetapeaks` output element list** — it omits `SummarySpectra` (a sixth return element) and does not name the sub-elements of `Untargeted` (`UntargetedIntensity`, `UntargetedCorrespondence`).
- [ ] **Correct the `preProcess` description** — README says it calls `Cardinal::smoothSpectra`, but the actual code calls `Cardinal::smooth`.
- [ ] **Fix the `assignMetapeaks` example argument order** — named args `refList = panel, pre = pre` are in the reverse order of the function signature `assignMetapeaks(x, pre, refList, ...)`.

### Add missing sections

- [ ] **Add a `## Citation` section** — point users to `citation("gutenTAG")` and acknowledge Cardinal (Bemis et al. 2023).
- [ ] **Add a `## License` section** (MIT, with link to LICENSE file).
- [ ] **Add Bioconductor installation instructions** (`BiocManager::install()`) as the primary route, with GitHub as the dev-version fallback.
- [ ] **Add system/dependency note** — Bioconductor dependencies cannot be installed via `install.packages()`; direct users to `BiocManager`.
- [ ] **Extend the Functionality section** to cover the QC and analysis steps promised in the introduction: `computeGearysC()`, `computeSNR()`, `computeVariogram()`, `louvainCluster()`, and the conversion functions (`asCytoImageList()`, `asSpatialExperiment()`, `asAnnData()`).
- [ ] **Add status badges** (R-CMD-check, Bioconductor build status, license).
- [ ] **Add `## Contributing` section** or link to `CONTRIBUTING.md` / issue tracker.

### Prose and formatting fixes

- [ ] **Fix typos in the opening paragraph**: "data handing" → "data handling"; "data is processing" → "data is processed".
- [ ] **Rewrite the opening paragraph** — it is a copy-paste of `DESCRIPTION`; the README should be a more readable, user-facing expansion.
- [ ] **Fix heading hierarchy** — numbered steps use `####` (H4) under a `##` (H2) parent; they should be `###` (H3).
- [ ] **Standardise code fence tags** — some use ` ```r ` (trailing space), others ` ```r `; make them uniform.
- [ ] **Fix the incomplete sentence** in the metapeak description: "contain information the same molecular species" should read "correspond to the same molecular species".
- [ ] **Remove the dangling `---` rule** at the end of the file — the document ends abruptly with a thematic break and no closing section.
- [ ] **Fix the `.ibd file` note** — currently says "not required to read in", which implies it's optional. Reword to clarify it is required for spectral data but does not need to be passed explicitly (Cardinal locates it automatically from the `.imzML` path).
- [ ] **Resolve `metapeakGeneration` vs `generateMetapeaks`** — the introduction and some roxygen blocks reference a non-existent function `metapeakGeneration`; the actual export is `generateMetapeaks`. Fix all occurrences in README prose.
- [ ] **Expand the `(and cytoviewer)` aside** — add a sentence explaining what cytoviewer is and how it relates to the package. (Leave wording to user.)
