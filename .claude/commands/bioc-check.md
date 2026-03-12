---
description: Full Bioconductor package health check — R CMD BUILD/INSTALL/CHECK, roxygen, dep-audit, and test coverage
allowed-tools: [Bash, Read, Glob, Grep]
---

# gutenTAG Bioconductor Health Check

Run a full package health check. Work through all phases below in order, printing a clearly labelled status line for each phase. Collect all warnings and errors and print a consolidated summary at the end.

## Phase 0 — Context

Run the following and record the values for use throughout:

```bash
git branch --show-current        # current branch → $BRANCH
Rscript -e "desc::desc_get('Package', file='DESCRIPTION')"   # package name
Rscript -e "desc::desc_get('Version', file='DESCRIPTION')"   # package version
```

If `$BRANCH` is `main` or `master`, this is a **release check** — apply stricter flags.
If `$BRANCH` is anything else, this is a **feature-branch check** — lighter flags, but still report all issues.

Print a header block:
```
=============================================================
  gutenTAG Bioconductor Health Check
  Branch : <branch>   Mode: <release|feature-branch>
  Package: <package>  Version: <version>
  Date   : <date>
=============================================================
```

## Phase 1 — roxygen2 / NAMESPACE sync

```r
devtools::document()
```

Check whether `NAMESPACE` or any `.Rd` files changed after running `devtools::document()`. If they did, report which files changed — this means the committed NAMESPACE/docs are out of sync with the roxygen tags.

## Phase 2 — Dependency audit

Scan all files under `R/` for `::` usage to find all packages referenced at runtime:

```bash
grep -roh '[A-Za-z][A-Za-z0-9.]*::' R/ | sed 's/:://' | sort -u
```

Compare the result against the `Imports:` and `Suggests:` fields in `DESCRIPTION`. Report:
- Packages used via `::` but **absent** from DESCRIPTION (missing deps)
- Packages declared in DESCRIPTION `Imports` but **never** used via `::` or `@importFrom` (unused deps)

## Phase 3 — R CMD INSTALL

```bash
R CMD INSTALL --no-multiarch --with-keep.source .
```

Capture stdout and stderr. Report SUCCESS or FAILED with any error lines.

## Phase 4 — R CMD BUILD

```bash
R CMD BUILD --no-build-vignettes .
```

This produces a tarball named `<Package>_<Version>.tar.gz`. Record the tarball path for Phase 5.

On a feature branch, add `--no-build-vignettes` (already included above) to keep builds fast.
On `main`/`master`, also build vignettes: use `R CMD BUILD .` without `--no-build-vignettes`.

## Phase 5 — R CMD CHECK

On **`main`/`master`** (release check):
```bash
R CMD CHECK --as-cran <tarball>
```

On a **feature branch**:
```bash
R CMD CHECK --no-vignettes --no-manual <tarball>
```

Parse the output and categorise findings into:
- ERRORs
- WARNINGs
- NOTEs

## Phase 6 — Bioconductor-specific checks

```r
BiocCheck::BiocCheck(package = ".", new.pkg = FALSE, quit.with.status = FALSE)
```

If `BiocCheck` is not installed, note this and skip. Report all findings.

## Phase 7 — Test coverage

```r
covr::package_coverage(quiet = FALSE)
covr::zero_coverage(cov)   # list untested functions
```

If `covr` is not installed, fall back to:
```r
devtools::test()
```

Report: overall coverage %, and any functions with 0% coverage.

## Final Summary

Print a consolidated report:

```
=============================================================
  SUMMARY
=============================================================
  Phase 1 roxygen sync     : PASS / WARN (N files changed)
  Phase 2 dep audit        : PASS / WARN (N missing, N unused)
  Phase 3 R CMD INSTALL    : PASS / FAIL
  Phase 4 R CMD BUILD      : PASS / FAIL
  Phase 5 R CMD CHECK      : PASS / N errors / N warnings / N notes
  Phase 6 BiocCheck        : PASS / N errors / N warnings / N notes
  Phase 7 Coverage         : XX% (N functions at 0%)
=============================================================
  Overall: PASS / FAIL
=============================================================
```

Then list all individual errors and warnings grouped by phase, with file and line references where available.

Clean up the tarball produced in Phase 4 after the check completes, unless the user passed `keep` as an argument (`$ARGUMENTS`).
