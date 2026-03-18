---

## Code Review Findings

### Priority 1 — Bugs that will break on clean install or produce silent wrong results


- [ ] **Determine minimum compatible Cardinal version**
  Currently pinned to `>= 3.6.2` in DESCRIPTION to avoid a version mismatch (Cardinal 3.x changed `summarizeFeatures` API: `stat=` → `FUN=`, dropping `"nnzero"`). Investigate the oldest Cardinal version that works correctly with gutenTAG and update the version pin accordingly.
  **Note:** Cardinal 3.6.2 can be installed on R 4.3.1 (without upgrading to R 4.4) by installing `matter >= 2.6.2` and `Cardinal` directly via `remotes::install_github()`. Standard `BiocManager::install()` on Bioc 3.17/3.18 is insufficient as it ships older versions of these packages.






### Priority 2 — Convention violations and R CMD CHECK failures

### Priority 3 — Performance and style issues

- [ ] **Add BiocParallel parallelisation to `generateMetapeaks` and `assignMetapeaks`**
  Plan saved at `~/.claude/plans/generic-bouncing-scone.md`. Parallelise the `for` loop over metapeaks in `estimateMetapeaks` (lines 37–76) and the `for` loop in `.generateFinalIntensityDF` (lines 136–141) using `BiocParallel::bplapply`. Add `BPPARAM = BiocParallel::bpparam()` to both exported functions and their internal helpers. Add validity checks and `@param BPPARAM` roxygen entries. Add serial-vs-parallel equality tests using `SnowParam(2)`.


### Priority 4 — Low-priority cleanup


